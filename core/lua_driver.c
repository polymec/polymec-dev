// Copyright (c) 2012-2019, Jeffrey N. Johnson
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

static bool _lua_driver_running = false;
static char* _script = NULL;

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
#define EOFMARK "<eof>"
#define marklen (sizeof(EOFMARK)/sizeof(char) - 1)

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

// This function writes an executable chunk to a string that can be
// sent across the wire. Adapted from 3rdparty/lua/src/lstrlib.c.
static int serialize_chunk(lua_State *L, const void *b, size_t size, void *B)
{
  luaL_addlstring((luaL_Buffer*)B, (const char *)b, size);
  return 0;
}

// This function broadcasts an executable chunk at the top of the stack
// from process 0 to all other processes.
static int broadcast_chunk(lua_State *L)
{
  ASSERT(_mpi_rank == 0);
  if (_mpi_nprocs > 1)
  {
    int strip = 0; // Don't bother stripping debug info.
    luaL_checktype(L, -1, LUA_TFUNCTION);

    luaL_Buffer buffer;
    luaL_buffinit(L, &buffer);
    if (lua_dump(L, serialize_chunk, &buffer, strip) != 0)
      return luaL_error(L, "Could not serialize compiled input on rank 0.");

    // Now broadcast the contents of the buffer to other ranks.
    int n = (int)buffer.n;
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(buffer.b, n, MPI_CHAR, 0, MPI_COMM_WORLD);

    // Consolidate the stack and pop the buffer off the top.
    luaL_pushresult(&buffer);
    lua_pop(L, 1);
  }
  return LUA_OK;
}

static void broadcast_error(lua_State *L, int err)
{
  ASSERT(_mpi_rank == 0);
  if (_mpi_nprocs > 1)
  {
    // Broadcast a -1 to signify that an error occurred.
    int n = -1;
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Now broadcast the error code.
    MPI_Bcast(&err, 1, MPI_INT, 0, MPI_COMM_WORLD);
  }
}

static int receive_chunk(lua_State *L)
{
  ASSERT(_mpi_rank != 0);
  ASSERT(_mpi_nprocs > 1);

  // Receive the broadcast from rank 0.
  int n;
  MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (n >= 0) // no error!
  {
    char* buffer = polymec_malloc(sizeof(char) * (n+1));
    MPI_Bcast(buffer, n, MPI_CHAR, 0, MPI_COMM_WORLD);
    buffer[n] = '\0';
    // Compile the buffer into an executable chunk.
    int status = luaL_loadbuffer(L, (const char*)buffer, (size_t)n, "input");
    polymec_free(buffer);
    return status;
  }
  else
  {
    // Get the error code and return it.
    int status;
    MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD);
    return status;
  }
}

extern int register_thirdparty_module(lua_State* L);

// This function controls the main loop, and is called in Lua's protected mode.
static int pmain(lua_State* L)
{
  // We're up and running!
  _lua_driver_running = true;

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
  register_thirdparty_module(L);
  if (register_types_and_modules != NULL)
    register_types_and_modules(L);

  // Run the input file, if we have one.
  if (filename != NULL)
  {
    // Load the file on rank 0, parse it into an executable chunk, and
    // broadcast this chunk to other ranks.
    int status;
    if (_mpi_rank == 0)
    {
      status = luaL_loadfile(L, (const char*)filename);
      if (status == LUA_OK)
        status = broadcast_chunk(L);
      else
        broadcast_error(L, status);
    }
    else
      status = receive_chunk(L);

    // Execute the chunk on every rank.
    _script = filename;
    if (status == LUA_OK)
      status = lua_pcall(L, 0, LUA_MULTRET, 0);
    if (status != LUA_OK)
    {
      report_error(L, status);
      if (!interactive)
      {
        lua_pushboolean(L, false);
        return 1;
      }
    }
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
  _script = NULL;
  _lua_driver_running = false;
  lua_pushboolean(L, true);
  return 1;
}

extern lua_State* polymec_lua_State(void);

extern noreturn void exit(int status);

static noreturn void usage(int argc, char** argv)
{
  // Fire up MPI real quick and fetch our rank.
  int prev_initialized, rank;
  MPI_Initialized(&prev_initialized);
  if (!prev_initialized)
    MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Explain usage only on rank 0.
  if (rank == 0)
  {
    printf("%s: usage:\n", argv[0]);
    printf("%s [flags] <input> [option1=val1 [option2=val2 [...]]]\n\n", argv[0]);
    printf("Flags are:\n");
    printf(" -i              Starts an interactive session (valid for 1 process only).\n");
    printf(" --help, -h      Displays help.\n\n");
    printf("Options are:\n");
    printf(" pause=N         Pauses for N seconds before executing, showing PID(s).\n");
    printf(" log=LEVEL       Sets the level of detail at which messages are logged.\n");
    printf("                 Case-insensitive levels in ascending verbosity are:\n");
    printf("                 off,urgent,info,detail,debug\n");
    printf("                 (Default is info)\n");
    printf(" log_mode=MODE   Specifies how messages are logged.\n");
    printf("                 Case-insensitive modes are:\n");
    printf("                 single <-- Logs messages to MPI rank 0\n");
    printf("                 all    <-- Logs messages to all MPI ranks\n");
    printf("                 P      <-- Logs messages to MPI rank P\n");
    printf(" log_file=PREFIX Logs messages to PREFIX.log.\n");
    printf("                 If log_mode=all, logs messages for MPI\n");
    printf("                 rank P in PREFIX.P.log.\n");
    printf("                 Existing log files are overwritten.\n");
    printf(" mpi_errors=VAL  Sets up an error handler for MPI:\n");
    printf("                 fatal  <-- MPI errors are fatal\n");
    printf("                 return <-- MPI errors return error codes\n");
    printf(" num_threads=N   Sets number of OpenMP threads to use.\n");
    printf(" timers=VAL      Enables or disables timers.\n");
    printf("                 Case-insensitive values are:\n");
    printf("                 1,true,yes,on     <-- enable\n");
    printf("                 (everything else) <-- disable\n");
    printf(" timer_file=PATH Specifies the file for the timer report if\n");
    printf("                 timers=1. Default: timer_report.txt\n");
    printf(" dl_paths=PATH   Sets path(s) to search for dynamically loaded libraries.\n");
    printf("                 PATH is a colon-delimited list of directories.\n\n");
    printf("You can specify other options as well. All options are made available\n");
    printf("in the options table, accessible in your input.\n");
  }

  // Shut down if we have to.
  if (!prev_initialized)
    MPI_Finalize();
  exit(0);
}

static void lua_error_handler(const char* message)
{
  lua_State* L = polymec_lua_State();
  luaL_error(L, "%s", message);
}

typedef void (*sighandler_t)(int);
static void handle_sigint_or_sigterm(int sig)
{
  // Exit.
  exit(sig);
}

int lua_driver(int argc,
               char** argv,
               int (*register_types_and_modules)(lua_State* L))
{
  // Look for "-h" or "--help" in the command line args and print usage
  // info if we find it.
  for (int i = 1; i < argc; ++i)
  {
    if ((strcmp(argv[i], "-h") == 0) || (strcmp(argv[i], "--help") == 0))
      usage(argc, argv);
  }

  // Start everything up.
  polymec_init(argc, argv);

  // Intercept SIGTERM and SIGINT and cause them to exit cleanly.
  sighandler_t def_sigterm_handler = signal(SIGTERM, handle_sigint_or_sigterm);
  sighandler_t def_sigint_handler = signal(SIGINT, handle_sigint_or_sigterm);

  // Take note of our rank and number of processes.
  MPI_Comm_rank(MPI_COMM_WORLD, &_mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &_mpi_nprocs);

  // Override the default error handler so that polymec_error reports to
  // Lua.
  polymec_set_error_handler(lua_error_handler);

  // Set the maximum history length.
  linenoiseHistorySetMaxLen(500);

  // Access our Lua state.
  lua_State* L = polymec_lua_State();

  // Pass the command line arguments into the main function, which runs in
  // protected mode.
  lua_pushcfunction(L, &pmain);
  lua_pushcfunction(L, register_types_and_modules);

  // Execute pmain, which parses the file and/or begins an interactive session.
  int status = lua_pcall(L, 1, 1, 0);
  int result = lua_toboolean(L, -1);

  // Reinstate the normal signal handlers.
  signal(SIGINT, def_sigint_handler);
  signal(SIGTERM, def_sigterm_handler);

  // Shut down and report any error(s).
  return (result && status == LUA_OK) ? EXIT_SUCCESS : EXIT_FAILURE;
}

bool lua_driver_running(void)
{
  return _lua_driver_running;
}

const char* lua_driver_script(void)
{
  return _script;
}

