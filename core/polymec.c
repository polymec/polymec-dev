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

#include <stdlib.h>
#include <time.h>
#include <stdarg.h>
#include <gc/gc.h>
#include "polytope_c.h"
#include "core/polymec.h"
#include "core/polymec_version.h"
#include "core/options.h"
#include "core/model.h"

#ifdef Linux
#include <fenv.h>
// FE_INEXACT    inexact result
// FE_DIVBYZERO  division by zero
// FE_UNDERFLOW  result not representable due to underflow
// FE_OVERFLOW   result not representable due to overflow
// FE_INVALID    invalid operation
#endif

#ifdef APPLE
#include <xmmintrin.h>
// _MM_MASK_DIV_INEXACT     inexact result
// _MM_MASK_DIV_ZERO        division by zero
// _MM_MASK_DIV_OVERFLOW    result not reproducible due to overflow
// _MM_MASK_DIV_UNDERFLOW   result not reproducible due to underflow
// _MM_MASK_INVALID         invalid operation
// _MM_MASK_DENORM          denormalized number
#endif

// Error handler.
static polymec_error_handler_function error_handler = NULL;

// Functions to call on exit.
typedef void (*at_exit_func)();
static at_exit_func _atexit_funcs[32];
static int _num_atexit_funcs = 0;

static void shutdown()
{
  for (int i = 0; i < _num_atexit_funcs; ++i)
    _atexit_funcs[i]();

  MPI_Finalize();
#ifndef NDEBUG
  polymec_disable_fpe_exceptions();
#endif
}

void polymec_init(int argc, char** argv)
{
#ifndef NDEBUG
  polymec_enable_fpe_exceptions();
#endif

  // Start up MPI.
  MPI_Init(&argc, &argv);

  // Start up the garbage collector.
  GC_INIT();

  // Register a shutdown function.
  polymec_atexit(&shutdown);
}

// Default error handler.
static int default_error_handler(const char* message)
{
  printf("Fatal error: %s\n", message);
//#if USE_MPI
//  MPI_Abort(MPI_COMM_WORLD, -1);
//#else
//  abort();
//#endif
  exit(-1);
  return -1; // Not reached.
}

int 
polymec_error(const char* message, ...)
{
  // Set the default error handler if no handler is set.
  if (error_handler == NULL)
    error_handler = &default_error_handler;

  // Extract the variadic arguments and splat them into a string.
  char err[1024];
  va_list argp;
  va_start(argp, message);
  vsnprintf(err, 1024, message, argp);
  va_end(argp);

  // Call the handler.
  return error_handler(err);
}

static void polymec_polytope_error_handler(const char* error, int status)
{
  polymec_error(error);
}

void 
polymec_set_error_handler(polymec_error_handler_function handler)
{
  error_handler = handler;
  polytope_set_error_handler(polymec_polytope_error_handler);
}

void 
polymec_warn(const char* message, ...)
{
  // Extract the variadic arguments and splat them into a string.
  va_list argp;
  va_start(argp, message);
  vfprintf(stderr, message, argp);
  va_end(argp);
}

void 
polymec_enable_fpe_exceptions()
{
#ifndef NDEBUG

#ifdef Linux
  feclearexcept(FE_ALL_EXCEPT);
  int flags = FE_DIVBYZERO |
              FE_INVALID   |
//              FE_UNDERFLOW |
              FE_OVERFLOW;

  feenableexcept(flags);
#endif

#ifdef APPLE
  unsigned int mask = _MM_MASK_INVALID | _MM_MASK_DIV_ZERO | _MM_MASK_OVERFLOW;
  _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() &~ mask);
#endif

#endif
}

void 
polymec_disable_fpe_exceptions()
{
#ifndef NDEBUG

#ifdef Linux
  fedisableexcept(fegetexcept());
#endif

#ifdef APPLE
  unsigned int mask = _MM_MASK_INVALID | _MM_MASK_DIV_ZERO | _MM_MASK_OVERFLOW;
  _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & mask);
#endif

#endif
}

void polymec_not_implemented(const char* component)
{
  char err[1024];
  snprintf(err, 1024, "polymec: not implemented: %s\n", component);
  default_error_handler(err);
}

void polymec_atexit(void (*func)()) 
{
  ASSERT(_num_atexit_funcs <= 32);
  _atexit_funcs[_num_atexit_funcs++] = func;
}

void polymec_version_fprintf(const char* exe_name, FILE* stream)
{
  fprintf(stream, "%s v%s\n", exe_name, POLYMEC_VERSION);
}

void polymec_provenance_fprintf(int argc, char** argv, FILE* stream)
{
  fprintf(stream, "=======================================================================\n");
  fprintf(stream, "                                Provenance:\n");
  fprintf(stream, "=======================================================================\n");
  fprintf(stream, "Version: %s\n", POLYMEC_VERSION);
  fprintf(stream, "Invoked with: ");
  for (int i = 0; i < argc; ++i)
    fprintf(stream, "%s ", argv[i]);
  fprintf(stream, "\n");
  time_t raw_time;
  time(&raw_time);
  fprintf(stream, "Invoked on: %s\n", ctime(&raw_time));

  if (strlen(POLYMEC_GIT_DIFF) > 0)
  {
    fprintf(stream, "=======================================================================\n");
    fprintf(stream, "Modifications to revision:\n");
    fprintf(stream, "%s\n\n", POLYMEC_GIT_DIFF);
  }

  // If we received an input script, write out its contents.
  options_t* options = options_parse(argc, argv);
  char* input = options_input(options);
  if (input == NULL)
  {
    // It's possible that the 1st argument is actually the input file.
    char* command = options_command(options);
    FILE* fp = fopen(command, "r");
    if (fp != NULL)
    {
      input = command;
      fclose(fp);
    }
  }
  if (input != NULL)
  {
    FILE* fp = fopen(input, "r");
    if (fp == NULL)
      fprintf(stream, "Invalid input specified.");
    else
    {
      char buff[sizeof(char)*1024];
      fseek(fp, 0L, SEEK_END);
      int end = ftell(fp);
      rewind(fp);
      static const int input_len_limit = 10 * 1024 * 1024; // 10kb input limit.
      bool truncated = false;
      if (input_len_limit < end)
      {
        truncated = true;
        end = input_len_limit; 
      }
      fprintf(stream, "=======================================================================\n");
      fprintf(stream, "Contents of input script:\n");
      int offset = 0;
      while (offset < end)
      {
        fread(buff, 1, MIN(1000, end-offset), fp);
        if (end-offset < 1000)
          buff[end-offset] = '\0';
        fprintf(stream, "%s", buff);
        offset += MIN(1000, end-offset);
      }
      if (truncated)
        fprintf(stream, "\n<<< truncated >>>\n");
      fclose(fp);
    }
    fprintf(stream, "\n");
  }
  options = NULL;

  fprintf(stream, "=======================================================================\n\n");
}

