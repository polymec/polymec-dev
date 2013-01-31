#include <stdlib.h>
#include <stdarg.h>
#include <gc/gc.h>
#include "core/polymec.h"
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

#ifdef __cplusplus
extern "C" {
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
  printf("Error: %s\n", message);
  printf("encountered in polymec. Exiting with status -1\n");
#if USE_MPI
  MPI_Abort(MPI_COMM_WORLD, -1);
  return -1; // Not reached.
#else
  abort();
#endif
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

void 
polymec_set_error_handler(polymec_error_handler_function handler)
{
  error_handler = handler;
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
  _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID
                                                  & ~_MM_MASK_DIV_ZERO
                                                  & ~_MM_MASK_OVERFLOW);
#endif
#endif
}

void 
polymec_disable_fpe_exceptions()
{
#ifdef Linux
#ifndef NDEBUG
  fedisableexcept(fegetexcept());
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

#ifdef __cplusplus
}
#endif
