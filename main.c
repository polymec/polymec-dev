#include "core/arbi.h"
#include "core/options.h"
#include "core/simulation.h"
#include "core/soln_vector.h"

#ifdef __cplusplus
extern "C" {
#endif

static void usage()
{
  fprintf(stderr, "arbi: usage:\n");
  fprintf(stderr, "arbi [options] <input_file>\n");
  exit(-1);
}

int main(int argc, char** argv)
{
  // Start everything up.
  arbi_init(argc, argv);

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

#ifdef __cplusplus
}
#endif
