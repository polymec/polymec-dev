#include "core/polymec.h"
#include "core/options.h"
#include "core/interpreter.h"
#include "geometry/interpreter_register_geometry_functions.h"

#ifdef __cplusplus
extern "C" {
#endif

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
extern int write_tough_mesh(lua_State* lua);
extern int write_vtk_plot(lua_State* lua);

static void interpreter_register_mesher_functions(interpreter_t* interpreter)
{
  interpreter_register_function(interpreter, "write_tough_mesh", write_tough_mesh);
  interpreter_register_function(interpreter, "write_vtk_plot", write_vtk_plot);
}

int main(int argc, char** argv)
{
  // Start everything up.
  polymec_init(argc, argv);

  // Set the log level.
  set_log_level(LOG_INFO);

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

#ifdef __cplusplus
}
#endif
