#include <stdlib.h>
#include "core/arbi.h"
#include "core/options.h"
#include "core/model.h"
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

  // Extract the name of the desired model.
  char* model_name = options_model(opts);

  // Attempt to construct the model.
  model_t* model = arbi_model(model_name, opts);

  // Have we been asked for help with the model?
  if (options_help(opts))
  {
    model_usage(model, stdout);
    exit(-1);
  }

  // Have we been asked to run a benchmark?
  char* benchmark = options_benchmark(opts);
  if (benchmark != NULL)
  {
    model_run_benchmark(model, benchmark);
    exit(0);
  }

  // Create a simulation in which to execute the model.
  simulation_t* sim = simulation_new(model, opts);

  // Initialize the simulation.
  double t = simulation_start_time(sim);
  simulation_init(sim, model, t);

  // Run the simulation.
  int num_steps = simulation_num_steps(sim);
  double max_time = simulation_max_time(sim);
  int step = 0;
  while ((t < max_time) && (step < num_steps))
  {
    simulation_invoke_callbacks(sim, model, t, step);
    simulation_step(sim, model, &t);
    ++step;
  }

  // Clean up.
  model_free(model);
  simulation_free(sim);
  options_free(opts);

  // That's it.
  return 0;
}

#ifdef __cplusplus
}
#endif
