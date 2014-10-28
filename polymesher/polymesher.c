// Copyright (c) 2012-2014, Jeffrey N. Johnson
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this 
// list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice, 
// this list of conditions and the following disclaimer in the documentation 
// and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include <strings.h>
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
  fprintf(stream, "\nType 'polymesher help' for documentation.\n");
}

static void mesher_help(interpreter_t* interp, const char* topic, FILE* stream)
{
  // If no argument was given, just print the polymesher's basic 
  // documentation.
  if (topic == NULL)
  {
    fprintf(stream, "polymesher: A polyhedral mesh generator.\n\n");
    fprintf(stream, "polymesher executes Lua scripts that create and manipulate polyhedral meshes\n");
    fprintf(stream, "using a variety of functions and objects.\n");
    fprintf(stream, "\nUse 'polymesher help list' to list available functions, and\n");
    fprintf(stream, "'polymesher help list <function>' for documentation on a given function.\n");
  }
  else
  {
    // Attempt to dig up the documentation for the given registered function.
    interpreter_help(interp, topic, stream);
  }
}

// Lua stuff.
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"

// Interpreter functions.
extern int write_gnuplot_points(lua_State* lua);

static void interpreter_register_mesher_functions(interpreter_t* interpreter)
{
  interpreter_register_function(interpreter, "write_gnuplot_points", write_gnuplot_points, NULL);
}

int main(int argc, char** argv)
{
  // Start everything up.
  polymec_init(argc, argv);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  options_t* opts = NULL;
  char* input = NULL;
  interpreter_t* interp = NULL;
  int leave = 0;
  if (rank == 0)
  {
    // Get the parsed command line options.
    opts = options_argv();

    // Extract the input file and arguments. 
    input = options_argument(opts, 1);
    if (input == NULL)
    {
      mesher_usage(stderr);
      leave = 1;
      goto leaving;
    }

    // Full provenance, or no?
    char* provenance_str = options_value(opts, "provenance");
    bool provenance = ((provenance_str != NULL) && !strcmp(provenance_str, "1"));

    if ((strcmp(input, "help") != 0) && (rank == 0))
    {
      // Check to see whether the given file exists.
      FILE* fp = fopen(input, "r");
      if (fp == NULL)
      {
        fprintf(stderr, "polymesher: Input file not found: %s\n", input);
        leave = 1;
        goto leaving;
      }
      else
        fclose(fp);
    }

    // Set the log level.
    log_level_t log_lev = LOG_DETAIL;
    char* logging = options_value(opts, "logging");
    if (logging != NULL)
    {
      if (!strcasecmp(logging, "debug"))
        log_lev = LOG_DEBUG;
      else if (!strcasecmp(logging, "detail"))
        log_lev = LOG_DETAIL;
      else if (!strcasecmp(logging, "info"))
        log_lev = LOG_INFO;
      else if (!strcasecmp(logging, "urgent"))
        log_lev = LOG_URGENT;
      else if (!strcasecmp(logging, "off"))
        log_lev = LOG_NONE;
    }
    set_log_level(log_lev);
    FILE* log_str = log_stream(log_lev);

    // Print a version identifier.

    // If we're providing full provenance, do so here.
    if (provenance)
      polymec_provenance_fprintf(log_str);
    else
      polymec_version_fprintf("polymesher", log_str);

    // Set up an interpreter for parsing the input file.
    interp = interpreter_new(NULL);
    interpreter_register_geometry_functions(interp);
    interpreter_register_mesher_functions(interp);

    // If we were asked for help, service the request here.
    if (!strcmp(input, "help"))
    {
      char* topic = options_argument(opts, 2);
      mesher_help(interp, topic, stderr);
      leave = 1;
      goto leaving;
    }
  }
  else
  {
    opts = options_argv();
    if (opts != NULL)
      input = options_argument(opts, 1);
    if (input != NULL)
    {
      interp = interpreter_new(NULL);
      interpreter_register_geometry_functions(interp);
      interpreter_register_mesher_functions(interp);
    }
  }

leaving:
  // Did something go wrong?
  MPI_Bcast(&leave, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if (leave)
    exit(0);

  // Parse it!
  interpreter_parse_file(interp, input);

  // Clean up.
  interpreter_free(interp);

  return 0;
}

