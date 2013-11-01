// Copyright 2012-2013 Jeffrey Johnson.
// 
// This file is part of Polymec, and is licensed under the Apache License, 
// Version 2.0 (the "License"); you may not use this file except in 
// compliance with the License. You may may find the text of the license in 
// the LICENSE file at the top-level source directory, or obtain a copy of 
// it at
// 
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "core/polymec.h"
#include "core/options.h"
#include "core/interpreter.h"
#include "geometry/interpreter_register_geometry_functions.h"

static void mesher_usage(FILE* stream)
{
  polymec_version_fprintf("polymesher", stream);
  fprintf(stream, "usage: polymesher [file] [options]\n\n");
  fprintf(stream, "Here, [file] is a file specifying instructions for generating a mesh.\n");
  fprintf(stream, "Options are:\n");
  fprintf(stream, "  provenance={*0*,1} - provides full provenance information (w/ diffs)\n");
  exit(-1);
}

// Lua stuff.
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"

// Interpreter functions.
extern int mesh_nodes(lua_State* lua);
extern int write_tough_mesh(lua_State* lua);
//extern int write_t2v_mesh(lua_State* lua);
extern int write_gnuplot_points(lua_State* lua);
//extern int read_meshvoro_mesh(lua_State* lua);

static void interpreter_register_mesher_functions(interpreter_t* interpreter)
{
  interpreter_register_function(interpreter, "mesh_nodes", mesh_nodes);
//  interpreter_register_function(interpreter, "read_meshvoro_mesh", read_meshvoro_mesh);
  interpreter_register_function(interpreter, "write_gnuplot_points", write_gnuplot_points);
//  interpreter_register_function(interpreter, "write_t2v_mesh", write_t2v_mesh);
  interpreter_register_function(interpreter, "write_tough_mesh", write_tough_mesh);
//  interpreter_register_function(interpreter, "read_meshvoro_mesh", read_meshvoro_mesh);
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

  // Full provenance, or no?
  char* provenance_str = options_value(opts, "provenance");
  bool provenance = ((provenance_str != NULL) && !strcmp(provenance_str, "1"));

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

  // Print a version identifier.
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0)
  {
    // If we're providing full provenance, do so here.
    if (provenance)
      polymec_provenance_fprintf(argc, argv, stdout);
    else
      polymec_version_fprintf("polymesher", stdout);
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

