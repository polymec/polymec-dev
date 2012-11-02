#include <string.h>
#include <stdlib.h>
#include "core/arbi.h"
#include "core/options.h"
#include "core/simulation.h"
#include "poisson/poisson_model.h"

#ifdef __cplusplus
extern "C" {
#endif

static const char* usage_str = 
"poisson: usage:\n"
"poisson [command] [args]\n\n"
"Here, [command]\n"
"is one of the following:\n\n"
"   run [filename]       -- Runs a simulation with input from the given file.\n"
"   benchmark [name]     -- Runs the given benchmark problem.\n"
"   help                 -- Prints information about the given model.\n\n";

static void usage()
{
  fprintf(stderr, "%s\n", usage_str);
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

  // Extract the command and arguments.
  char* command = options_command(opts);
  char* input = options_input(opts);
  if (!strcmp(command, "help"))
    usage();

  // Validate our inputs.
  if (command == NULL)
  {
    fprintf(stderr, "arbi: no command given! Usage:\n");
    fprintf(stderr, "arbi [command] [command args]\n");
    exit(-1);
  }
  static const char* valid_commands[] = {"run", "benchmark", "help", NULL};
  int c = 0;
  while (valid_commands[c] != NULL)
  {
    if (!strcmp(command, valid_commands[c]))
      break;
    ++c;
  }
  if (valid_commands[c] == NULL)
  {
    fprintf(stderr, "arbi: invalid command: '%s'\n", command);
    exit(-1);
  }

  // Attempt to construct the model.
  model_t* model = poisson_model_new(opts);
  ASSERT(model != NULL);

  // Have we been asked to run a benchmark?
  if (!strcmp(command, "benchmark"))
  {
    if (input == NULL)
    {
      fprintf(stderr, "poisson: No benchmark problem given! Usage:\n");
      fprintf(stderr, "poisson benchmark [problem]\n");
      exit(-1);
    }
    model_run_benchmark(model, input);
    exit(0);
  }

  // We are asked to run a simulation.
  ASSERT(!strcmp(command, "run"));
  if (input == NULL)
  {
    fprintf(stderr, "poisson: No input file given! Usage:\n");
    fprintf(stderr, "poisson run [input file]\n");
    exit(-1);
  }

  // Create a simulation in which to execute the model.
  simulation_t* sim = simulation_new(model, input, opts);

  // Initialize the simulation.
  simulation_init(sim);

  // Run the simulation.
  simulation_run(sim);

  // Clean up.
  simulation_free(sim);
  model_free(model);

  // That's it.
  return 0;
}

#ifdef __cplusplus
}
#endif
