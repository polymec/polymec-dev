#ifndef CLIDE_H
#define CLIDE_H

#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <mpi.h>

// Some macros.
#ifndef NDEBUG
#define ASSERT(x) \
  if (!(x)) \
  { \
    printf("Assertion %s failed\nat %s: %d\n", #x, __FILE__, __LINE__); \
    MPI_Abort(MPI_COMM_WORLD, -1); \
  }
#else
#define ASSERT(x) 
#endif

// This macro allows us to explicate that a function argument is not used
// in a C-friendly way. Some people point out that using this macro with 
// a "volatile" variable changes the structure of memory barriers. Perhaps 
// now you understand why clide isn't intended to be thread-safe.
#ifndef UNUSED_ARG
#define UNUSED_ARG(x) (void)(x)
#endif

// Here's a macro that fills an array with the given values.
#define FILL_ARRAY(array, start, num, val) \
for (int i = start; i < start + num; ++i) \
  array[i] = val

#define real_t double

// Error codes.
#define CLIDE_SUCCESS 0
#define CLIDE_FAILURE -1
#define CLIDE_NO_EFFECT 1

#ifdef __cplusplus
extern "C" {
#endif

typedef void (*clide_error_handler_function)(const char*);

// Issues an error with the given message. By default, an error issues a 
// message to stdout and exits the program with status -1.
void clide_error(const char* message);

// Sets the error handler for the clide library.
void clide_set_error_handler(clide_error_handler_function handler);

// Issues a warning to stderr.
void clide_warning(const char* message);

// This function enables floating point exceptions where available.
void clide_enable_fpe_exceptions();

//! This function disables floating point exceptions.
void clide_disable_fpe_exceptions();

#ifdef __cplusplus
}
#endif

#endif
