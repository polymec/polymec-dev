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

#include <unistd.h>
#include <time.h>
#include <stdarg.h>
#include <gc/gc.h>
#include "arena/proto.h"
#include "arena/pool.h"
#include "core/polymec.h"
#include "core/polymec_version.h"
#include "core/options.h"

// Standard C support for floating point environment.
#include <fenv.h>

#ifdef Linux
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

// This flag tells us whether polymec is already initialized.
static bool polymec_initialized = false;

// This function must be called to take advantage of Jonathan Shewchuk's 
// exact geometric predicates.
extern void exactinit();

// Command line arguments (used for provenance information).
int polymec_argc = 0;
char** polymec_argv = NULL;

// Invocation string and time.
char* polymec_invoc_str = NULL;
time_t polymec_invoc_time = 0;

// Error handler.
static polymec_error_handler_function error_handler = NULL;

// Functions to call on exit.
typedef void (*at_exit_func)();
static at_exit_func _atexit_funcs[32];
static int _num_atexit_funcs = 0;

static void shutdown()
{
  ASSERT(polymec_initialized);

  // Kill command line arguments.
  free(polymec_invoc_str);
  for (int i = 0; i < polymec_argc; ++i)
    free(polymec_argv[i]);
  free(polymec_argv);

  // Call shutdown functions.
  for (int i = 0; i < _num_atexit_funcs; ++i)
    _atexit_funcs[i]();

  MPI_Finalize();
#ifndef NDEBUG
  polymec_disable_fpe_exceptions();
#endif
}

void polymec_init(int argc, char** argv)
{
  if (!polymec_initialized)
  {
    // Jot down the invocation time.
    polymec_invoc_time = time(NULL);

#ifndef NDEBUG
    // By default, we enable floating point exceptions for debug builds.
    polymec_enable_fpe_exceptions();
#endif

    // Jot down command line args (use regular malloc).
    polymec_argc = argc;
    polymec_argv = malloc(sizeof(char*) * argc);
    for (int i = 0; i < argc; ++i)
      polymec_argv[i] = string_dup(argv[i]);

    // Construct the invocation string.
    int invoc_len = 2;
    for (int i = 0; i < polymec_argc; ++i)
      invoc_len += 1 + strlen(polymec_argv[i]);
    polymec_invoc_str = malloc(sizeof(char) * invoc_len);
    polymec_invoc_str[0] = '\0';
    for (int i = 0; i < polymec_argc; ++i)
    {
      char arg[strlen(polymec_argv[i])+2];
      if (i == polymec_argc-1)
        sprintf(arg, "%s", polymec_argv[i]); 
      else
        sprintf(arg, "%s ", polymec_argv[i]);
      strcat(polymec_invoc_str, arg);
    }

    // Start up MPI.
    MPI_Init(&argc, &argv);

    // Start up the garbage collector.
    GC_INIT();

    // Initialize variables for exact arithmetic.
    exactinit();

    // Register a shutdown function.
    atexit(&shutdown);

    // Okay! We're initialized.
    polymec_initialized = true;
  }
}

// Default error handler.
static noreturn void default_error_handler(const char* message)
{
  printf("Fatal error: %s\n", message);
//#if USE_MPI
//  MPI_Abort(MPI_COMM_WORLD, -1);
//#else
//  abort();
//#endif
  exit(-1);
}

#ifdef POLYMEC_HAVE_MPI
// Here are the error handlers for the Scotch partitioning library.
void SCOTCH_errorPrint(const char* const errstr, ...)
{
  va_list argp;
  va_start(argp, errstr);
  polymec_error(errstr, argp);
  va_end(argp);
}

void SCOTCH_errorPrintW(const char* const errstr, ...)
{
  va_list argp;
  va_start(argp, errstr);
  polymec_warn(errstr, argp);
  va_end(argp);
}
#endif

void polymec_abort(const char* message, ...)
{
  // Extract the variadic arguments and splat them into a string.
  char err[1024];
  va_list argp;
  va_start(argp, message);
  vsnprintf(err, 1024, message, argp);
  va_end(argp);

  // Abort.
  fprintf(stderr, "%s\n", err);
#if USE_MPI
  MPI_Abort(MPI_COMM_WORLD, -1); 
#else
  abort(); 
#endif
}

void polymec_error(const char* message, ...)
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
  error_handler(err);

  // Make sure we don't return.
  exit(-1);
}

void polymec_set_error_handler(polymec_error_handler_function handler)
{
  error_handler = handler;
}

void polymec_warn(const char* message, ...)
{
  // Extract the variadic arguments and splat them into a string.
  va_list argp;
  va_start(argp, message);
  vfprintf(stderr, message, argp);
  va_end(argp);
}

void polymec_enable_fpe_exceptions()
{
  feclearexcept(FE_ALL_EXCEPT);
#ifdef Linux
  int flags = FE_DIVBYZERO |
              FE_INVALID   |
//              FE_UNDERFLOW |
              FE_OVERFLOW;

  feenableexcept(flags);
#endif

#ifdef APPLE
  unsigned int mask = _MM_MASK_INVALID | _MM_MASK_DIV_ZERO | _MM_MASK_OVERFLOW | _MM_MASK_DENORM;
  _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~mask);
#endif
}

void polymec_disable_fpe_exceptions()
{
#ifdef Linux
  fedisableexcept(fegetexcept());
#endif

#ifdef APPLE
  unsigned int mask = _MM_MASK_INVALID | _MM_MASK_DIV_ZERO | _MM_MASK_OVERFLOW | _MM_MASK_DENORM;
  _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & mask);
#endif
}


// The following suspend/restore mechanism uses standard C99 floating point 
// kung fu.
static fenv_t fpe_env;

void polymec_suspend_fpe_exceptions()
{
  // Hold exceptions till further notice.
  feholdexcept(&fpe_env);
}

void polymec_restore_fpe_exceptions()
{
  // Clear all exception flags.
  feclearexcept(FE_ALL_EXCEPT);

  // Now restore the previous floating point environment.
  fesetenv(&fpe_env);
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
  if (stream == NULL) return;
  fprintf(stream, "%s v%s\n", exe_name, POLYMEC_VERSION);
}

void polymec_provenance_fprintf(FILE* stream)
{
  ASSERT(polymec_initialized);
  if (stream == NULL) return;

  fprintf(stream, "=======================================================================\n");
  fprintf(stream, "                                Provenance:\n");
  fprintf(stream, "=======================================================================\n");
  fprintf(stream, "Version: %s\n", POLYMEC_VERSION);
  fprintf(stream, "Invoked with: %s\n", polymec_invoc_str);
  fprintf(stream, "Invoked on: %s\n", ctime(&polymec_invoc_time));

  if (POLYMEC_NUM_GIT_DIFFS > 0)
  {
    fprintf(stream, "=======================================================================\n");
    fprintf(stream, "Modifications to revision:\n");
    for (int i = 0; i < POLYMEC_NUM_GIT_DIFFS; ++i)
      fprintf(stream, "%s", POLYMEC_GIT_DIFFS[i]);
    fprintf(stream, "\n\n");
  }

  // If we received an input script, write out its contents.
  options_t* options = options_parse(polymec_argc, polymec_argv);
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
      char buff[1024];
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

const char* polymec_invocation()
{
  ASSERT(polymec_initialized);
  return (const char*)polymec_invoc_str;
}

time_t polymec_invocation_time()
{
  return polymec_invoc_time;
}

int polymec_num_cores()
{
#ifdef LINUX
  return sysconf(_SC_NPROCESSORS_ONLN);
#elif defined APPLE
  return sysconf(_SC_NPROCESSORS_ONLN);
#endif
}

