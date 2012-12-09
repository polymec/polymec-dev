#ifndef POLYMEC_H
#define POLYMEC_H

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <mpi.h>
#include "core/loggers.h"

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
// now you understand why polymec isn't intended to be thread-safe.
#ifndef UNUSED_ARG
#define UNUSED_ARG(x) (void)(x)
#endif

// Here's a macro that fills an array with the given values.
#define FILL_ARRAY(array, start, num, val) \
for (int i = start; i < start + num; ++i) \
  array[i] = val

#define real_t double

// This returns the minimum of a and b.
#ifndef MIN
#define MIN(a, b) ((a < b) ? a : b)
#endif

// This returns the maximum of a and b.
#ifndef MAX
#define MAX(a, b) ((a > b) ? a : b)
#endif

// This returns the sign (+1/-1) of the given quantity x.
#define SIGN(x) ((x > 0) ? 1 : -1)

// These macros use the given arena to allocate/free memory and fall back 
// to the heap if the first argument is NULL. Note that you still have to 
// #include <arena/arena.h> to use them.
#define ARENA_MALLOC(arena, size, alignment) \
  (arena != NULL) ? arena_malloc(arena, size, alignment) : malloc(size)

#define ARENA_REALLOC(arena, memory, size, alignment) \
  (arena != NULL) ? arena_realloc(arena, memory, size, alignment) : realloc(memory, size)

#define ARENA_FREE(arena, memory) \
  (arena != NULL) ? arena_free(arena, memory) : free(memory)

// Error codes.
#define POLYMEC_SUCCESS 0
#define POLYMEC_FAILURE -1
#define POLYMEC_NO_EFFECT 1

#ifdef __cplusplus
extern "C" {
#endif

typedef int (*polymec_error_handler_function)(const char*);

// Issues an error with the given message. By default, an error issues a 
// message to stdout and exits the program with status -1.
int polymec_error(const char* message, ...);

// Sets the error handler for the polymec library.
void polymec_set_error_handler(polymec_error_handler_function handler);

// Issues a warning to stderr.
void polymec_warn(const char* message, ...);

// This function enables floating point exceptions where available.
void polymec_enable_fpe_exceptions();

// This function disables floating point exceptions.
void polymec_disable_fpe_exceptions();

// Call this function to indicate that something is not implemented
// with a fatal exit.
void polymec_not_implemented(const char* component);

// Initialize polymec's services.
void polymec_init(int argc, char** argv);

// Register a function to be called upon exit.
void polymec_atexit(void (*func)());

#ifdef __cplusplus
}
#endif

#endif
