// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <signal.h>

#include "mpi.h"

#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"

#include "core/polymec.h"
#include "core/options.h"
#include "core/lua_driver.h"

// We use linenoise since readline has a complicated license, and libedit
// has a complicated dependency.
#include "linenoise.h"
#define lua_readline(L, b, p) ((void)L, ((b)=linenoise(p)) != NULL)
#define lua_saveline(L, line) ((void)L, linenoiseHistoryAdd(line))
#define lua_freeline(L, b)    ((void)L, free(b))

//------------------------------------------------------------------------
// The contents of this file are taken from a modified version of lua.c in 
// the Lua distribution.
//------------------------------------------------------------------------

static int _mpi_rank = -1;
static int _mpi_nprocs = -1;

// This function reports an error in parsing input.
static int report_error(lua_State *L, int status) 
{
  // Errors are only reported on MPI rank 0.
  if ((status != LUA_OK) && (_mpi_rank == 0))
  {
    const char *msg = lua_tostring(L, -1);
    lua_writestringerror("%s: ", polymec_executable_name());
    lua_writestringerror("%s\n", msg);
    lua_pop(L, 1); 
  }
  return status;
}

// This function handles error messages.
static int handle_message(lua_State *L) 
{
  const char *msg = lua_tostring(L, 1);
  if (msg == NULL) // not a string!
  {
    // Try to call __tostring on the message to render a string.
    if (luaL_callmeta(L, 1, "__tostring") && 
        (lua_type(L, -1) == LUA_TSTRING)) 
      return 1; 
    else
    {
      // We're stuck with a non-string, so barf out the type.
      msg = lua_pushfstring(L, "(error object is a %s value)",
                            luaL_typename(L, 1));
    }
  }

  // Return a standard traceback.
  luaL_traceback(L, L, msg, 1); 
  return 1;
}

// This function halts the interpreter when an interrupt signal is intercepted.
static void stop_interpreter(lua_State *L, lua_Debug *ar) 
{
  (void)ar; // unused arg

  // Reset the debug hook.
  lua_sethook(L, NULL, 0, 0);

  // Issue an error mesage.
  luaL_error(L, "Interrupted.");
}

// Global interpreter state for handling interrupts.
static lua_State* global_L = NULL;

// This function handles the interrupt signal.
static void handle_interrupt(int i) 
{
  signal(i, SIG_DFL); // If another SIGINT happens, terminate process.
  lua_sethook(global_L, stop_interpreter, LUA_MASKCALL | LUA_MASKRET | LUA_MASKCOUNT, 1);
}

// This function executes a chunk of Lua code.
static int execute_chunk(lua_State *L, int n_arg, int n_res) 
{
  // Put a message handler in place for tracebacks.
  int base = lua_gettop(L) - n_arg;
  lua_pushcfunction(L, handle_message);
  lua_insert(L, base);

  // Execute the code, setting a signal handler for the duration.
  global_L = L;  
  signal(SIGINT, handle_interrupt); 
  int status = lua_pcall(L, n_arg, n_res, base);
  signal(SIGINT, SIG_DFL);
  lua_remove(L, base);
  global_L = NULL;

  return status;
}

// mark in error messages for incomplete statements 
#define EOFMARK		"<eof>"
#define marklen		(sizeof(EOFMARK)/sizeof(char) - 1)

// This function checks whether status signals a syntax error and the 
// error message at the top of the stack ends with the above mark for
// incomplete statements.
static int incomplete(lua_State *L, int status) 
{
  if (status == LUA_ERRSYNTAX) 
  {
    size_t lmsg;
    const char *msg = lua_tolstring(L, -1, &lmsg);
    if ((lmsg >= marklen) && 
        (strcmp(msg + lmsg - marklen, EOFMARK) == 0)) 
    {
      lua_pop(L, 1);
      return 1;
    }
  }
  return 0;  
}

// This function prompts the user for a line of Lua code in interative mode.
static bool prompt_for_line(lua_State *L) 
{
  // Assemble a command prompt and stick it on the stack.
  const char *prog_name = polymec_executable_name();
  size_t pnl = strlen(prog_name);
  char prompt[pnl+3];
  strncpy(prompt, prog_name, sizeof(char)*pnl);
  prompt[pnl] = '>';
  prompt[pnl+1] = ' ';
  prompt[pnl+2] = '\0';
  lua_pushstring(L, prompt);

  // Read the line of input.
  char buffer[512];
  char *b = buffer;
  bool got_input = lua_readline(L, b, prompt);

  // If we didn't get anything, bug out (prompt removed by caller).
  if (!got_input) 
    return false;

  lua_pop(L, 1); // remove prompt from the stack.
  size_t l = strlen(b);

  // Remove any newline from the end of the line.
  if ((l > 0) && (b[l-1] == '\n'))
    b[--l] = '\0'; 

  lua_pushlstring(L, b, l);
  lua_freeline(L, b);
  return true;
}

// This function prepends a return command in front of an expression in 
// an attempt to compile it.
static int prepend_return(lua_State *L) 
{
  const char *line = lua_tostring(L, -1); // original line
  const char *mod_line = lua_pushfstring(L, "return %s;", line);
  int status = luaL_loadbuffer(L, mod_line, strlen(mod_line), "=stdin");
  if (status == LUA_OK) 
  {
    lua_remove(L, -2);  // remove modified line 
    if (line[0] != '\0') // not empty?
      lua_saveline(L, line);  // keep history 
  }
  else
    lua_pop(L, 2);  // pop result from 'luaL_loadbuffer' and modified line 
  return status;
}

// This function reads multiple lines until a complete Lua statement 
// is assembled.
static int compile_statement(lua_State *L) 
{
  while (true)
  {  
    // Repeat until a complete statement is assembled.
    size_t len;
    const char *line = lua_tolstring(L, 1, &len); 
    int status = luaL_loadbuffer(L, line, len, "=stdin"); // try it 
    if (!incomplete(L, status) || !prompt_for_line(L)) 
    {
      lua_saveline(L, line); // keep history 
      return status; // cannot or should not try to add continuation line 
    }
    lua_pushliteral(L, "\n");  // add newline... 
    lua_insert(L, -2);  // ...between the two lines 
    lua_concat(L, 3);  // join them 
  }
}

// This function reads a line of Lua code and trys to compile it as 
// (1) an expression, and (2) a statement (if (1) fails). 
static int read_line(lua_State *L) 
{
  lua_settop(L, 0);
  if (!prompt_for_line(L))
    return -1;  // no input 
  int status;
  if ((status = prepend_return(L)) != LUA_OK)  // 'return ...' didn't work.
    status = compile_statement(L); // try as statement, maybe with continuation lines 
  lua_remove(L, 1); // remove line from the stack 
  lua_assert(lua_gettop(L) == 1);
  return status;
}

// This function prints the result of an expression.
static void print_result(lua_State *L) 
{
  int n = lua_gettop(L);
  if (n > 0) 
  { 
    luaL_checkstack(L, LUA_MINSTACK, "too many results to print");
    lua_getglobal(L, "print");
    lua_insert(L, 1);
    if (lua_pcall(L, n, 0, 0) != LUA_OK)
      log_urgent("error calling 'print' (%s)", lua_tostring(L, -1));
  }
}

// This function handles interactive input.
static void interact(lua_State* L)
{
  int status;
  do
  {
    status = read_line(L);
    if (status == LUA_OK)
      status = execute_chunk(L, 0, LUA_MULTRET);
    if (status == LUA_OK) 
      print_result(L);
    else 
      report_error(L, status);
  }
  while (status != -1);
  
  // Clear the stack.
  lua_settop(L, 0); 
  lua_writeline();
}

// This function controls the main loop, and is called in Lua's protected mode.
static int pmain(lua_State* L)
{
  int (*register_types_and_modules)(lua_State* L) = lua_tocfunction(L, 1);

  // Get the parsed command line options.
  options_t* opts = options_argv();
  bool interactive = options_has_argument(opts, "-i");

  // Try to identify a file in the argument list.
  char* filename = NULL;
  size_t num_args = options_num_arguments(opts);
  for (size_t i = 1; i < num_args; ++i)
  {
    char* arg = options_argument(opts, i);

    // Disregard named values and flags.
    if (!string_contains(arg, "=") && (arg[0] != '-')) 
      filename = arg;
  }

  // Tell the standard Lua libraries to ignore environment variables.
  lua_pushboolean(L, true);
  lua_setfield(L, LUA_REGISTRYINDEX, "LUA_NOENV");

  // Open standard libraries, and register extra types and modules.
  luaL_openlibs(L);
  if (register_types_and_modules != NULL)
    register_types_and_modules(L);

  // Run the input file, if we have one.
  if (filename != NULL)
  {
    int status = luaL_loadfile(L, (const char*)filename);
    if (status == LUA_OK)
      status = lua_pcall(L, 0, LUA_MULTRET, 0);
    if (status != LUA_OK)
      report_error(L, status);
  }

  // If we're interactive, surrender control.
  if (interactive || (filename == NULL))
  {
    // We don't currently support interactive mode on more than one 
    // process.
    if (_mpi_nprocs > 1)
      polymec_error("Interactive mode is not supported in a parallel environment.");

    interact(L);
  }

  // We made it to the end without incident.
  lua_pushboolean(L, true);
  return 1;
}

int lua_driver(int argc,
               char** argv,
               int (*register_types_and_modules)(lua_State* L))
{
  // Start everything up.
  polymec_init(argc, argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &_mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &_mpi_nprocs);

  // Set the maximum history length.
  linenoiseHistorySetMaxLen(500);

  // Create a Lua state.
  lua_State* L = luaL_newstate();
  if (L == NULL)
    polymec_error("%s: cannot create Lua interpreter: not enough memory.", argv[0]);

  // Pass the command line arguments into the main function, which runs in 
  // protected mode.
  lua_pushcfunction(L, &pmain);
  lua_pushcfunction(L, register_types_and_modules);

  // Execute pmain, which parses the file and/or begins an interactive session.
  int status = lua_pcall(L, 1, 1, 0);
  int result = lua_toboolean(L, -1);

  // Shut down and report any error(s).
  lua_close(L);
  return (result && status == LUA_OK) ? EXIT_SUCCESS : EXIT_FAILURE;
}
