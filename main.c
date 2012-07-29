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
  model_t* model = NULL;
  if (model_exists(model_name))
  {
    model = model_new(model_name, opts);
  }
  else
  {
    char err[1024];
    snprintf(err, 1024, "Invalid model: '%s'", model_name);
    arbi_error(err);
  }

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
  simulation_init(sim, t);

  // Run the simulation.
  simulation_run(sim);

  // Clean up.
  simulation_free(sim);
  model_free(model);
  options_free(opts);

  // That's it.
  return 0;
}

#ifdef __cplusplus
}
#endif
