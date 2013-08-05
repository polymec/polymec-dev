#include "core/polymec.h"
#include "core/options.h"
#include "core/interpreter.h"
#include "geometry/interpreter_register_geometry_functions.h"

static void mesher_usage(FILE* stream)
{
  fprintf(stream, "polymesher: usage:\n");
  fprintf(stream, "polymesher [file] [args]\n\n");
  fprintf(stream, "Here, [file] is a file specifying instructions\n");
  fprintf(stream, "for generating a mesh.\n\n");
  exit(-1);
}

// Lua stuff.
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"

// Interpreter functions.
extern int mesh_nodes(lua_State* lua);
extern int write_tough_mesh(lua_State* lua);
extern int write_silo_plot(lua_State* lua);
extern int write_vtk_plot(lua_State* lua);
extern int read_meshvoro_mesh(lua_State* lua);

static void interpreter_register_mesher_functions(interpreter_t* interpreter)
{
  interpreter_register_function(interpreter, "mesh_nodes", mesh_nodes);
  interpreter_register_function(interpreter, "read_meshvoro_mesh", read_meshvoro_mesh);
  interpreter_register_function(interpreter, "write_silo_plot", write_silo_plot);
  interpreter_register_function(interpreter, "write_vtk_plot", write_vtk_plot);
  interpreter_register_function(interpreter, "read_meshvoro_mesh", read_meshvoro_mesh);
}

int main(int argc, char** argv)
{
  // Start everything up.
  polymec_init(argc, argv);

  // Parse options on the command line.
  options_t* opts = options_parse(argc, argv);
  if (opts == NULL)
    mesher_usage(stderr);

  // Extract the input file and arguments. Note that we use "command" here 
  // to get the input file, since it comes first.
  char* input = options_command(opts);
  if (!strcmp(input, "help") || (input == NULL))
    mesher_usage(stderr);

  // Check to see whether the given file exists.
  FILE* fp = fopen(input, "r");
  if (fp == NULL)
  {
    fprintf(stderr, "polymesher: Input file not found: %s\n", input);
    return -1;
  }
  fclose(fp);

  // Set the log level.
  set_log_level(LOG_DETAIL);
  char* logging = options_value(opts, "logging");
  if (logging != NULL)
  {
    if (!strcasecmp(logging, "debug"))
      set_log_level(LOG_DEBUG);
    else if (!strcasecmp(logging, "detail"))
      set_log_level(LOG_DETAIL);
    else if (!strcasecmp(logging, "info"))
      set_log_level(LOG_INFO);
    else if (!strcasecmp(logging, "urgent"))
      set_log_level(LOG_URGENT);
    else if (!strcasecmp(logging, "off"))
      set_log_level(LOG_NONE);
  }

  // Set up an interpreter for parsing the input file.
  interpreter_t* interp = interpreter_new(NULL);
  interpreter_register_geometry_functions(interp);
  interpreter_register_mesher_functions(interp);

  // Parse it!
  interpreter_parse_file(interp, input);

  // Clean up.
  interpreter_free(interp);

  return 0;
}

