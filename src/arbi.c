#include <stdlib.h>
#include <mpi.h>
#include "arbi.h"
#include "options.h"
#include "simulation.h"
#include "soln_vector.h"

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
static arbi_error_handler_function error_handler = NULL;

static void shutdown()
{
  MPI_Finalize();
}

static void usage()
{
  fprintf(stderr, "arbi: usage:\n");
  fprintf(stderr, "arbi [options] <input_file>\n");
  exit(-1);
}

int main(int argc, char** argv)
{
  // Start everything up.
  MPI_Init(&argc, &argv);

  // Register a shutdown function.
  atexit(&shutdown);

  // Parse options on the command line.
  options_t* opts = options_parse(argc, argv);
  if (opts == NULL)
    usage();

  // If we are given an input file, read it.
  simulation_t* sim = NULL;
  char* filename = options_file(opts);
  if (filename != NULL)
  {
    FILE* file = fopen(filename, "r");
    sim = simulation_from_file(file, opts);
    fclose(file);
    ASSERT(sim != NULL);
  }
  else
  {
    ASSERT(false);
  }
  options_free(opts);

  // Initialize the problem.
  soln_vector_t* solution = simulation_create_vector(sim);
  double t = simulation_start_time(sim);
  if (simulation_init(sim, solution, t) != ARBI_SUCCESS)
    exit(-1);

  // Run the simulation.
  int num_steps = simulation_num_steps(sim);
  double max_time = simulation_max_time(sim);
  int step = 0;
  while ((t < max_time) && (step < num_steps))
  {
    simulation_invoke_callbacks(sim, solution, t, step);
    if (simulation_step(sim, solution, &t) != ARBI_SUCCESS)
      exit(-1);
    ++step;
  }

  // Clean up.
  soln_vector_free(solution);
  simulation_free(sim);

  // That's it.
  return 0;
}

// Default error handler.
static void default_error_handler(const char* message)
{
  printf("Error: %s\n", message);
  printf("encountered in arbi. Exiting with status -1\n");
  MPI_Abort(MPI_COMM_WORLD, -1);
}

void 
arbi_error(const char* message)
{
  // Set the default error handler if no handler is set.
  if (error_handler == NULL)
    error_handler = &default_error_handler;

  // Call it.
  error_handler(message);
}

void 
arbi_set_error_handler(arbi_error_handler_function handler)
{
  error_handler = handler;
}

void 
arbi_warn(const char* message)
{
  fprintf(stderr, "%s\n", message);
}

void 
arbi_enable_fpe_exceptions()
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
arbi_disable_fpe_exceptions()
{
#ifdef Linux
#ifndef NDEBUG
  fedisableexcept(fegetexcept());
#endif
#endif
}

void arbi_not_implemented(const char* component)
{
  char err[1024];
  snprintf(err, 1024, "arbi: not implemented: %s\n", component);
  default_error_handler(err);
}

#ifdef __cplusplus
}
#endif
