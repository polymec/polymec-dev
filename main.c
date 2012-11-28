#include <string.h>
#include <stdlib.h>
#include "core/polymec.h"
#include "core/options.h"
#include "core/model.h"
#include "core/simulation.h"

#ifdef __cplusplus
extern "C" {
#endif

static const char* usage_str = 
"polymec: usage:\n"
"polymec [model] [command] [args] OR\n"
"polymec list\n\n"
"Here, [model] is a name that identifies a numerical model, and [command]\n"
"is one of the following:\n\n"
"   run [filename]       -- Runs a simulation with input from the given file.\n"
"   benchmark [name]     -- Runs the given benchmark problem.\n"
"   help                 -- Prints information about the given model.\n\n"
"If the list command is given, a list of available models is written to stdout.\n\n";

static void usage()
{
  fprintf(stderr, "%s\n", usage_str);
  exit(-1);
}

int main(int argc, char** argv)
{
  // Start everything up.
  polymec_init(argc, argv);

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
    fprintf(stderr, "polymec: no command given! Usage:\n");
    fprintf(stderr, "polymec [command] [command args]\n");
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
    fprintf(stderr, "polymec: invalid command: '%s'\n", command);
    exit(-1);
  }

#if 0
  // Attempt to construct the model.
  model_t* model = model_new(model_name, opts);
  ASSERT(model != NULL);

  // Have we been asked for help with the model?
  if (!strcmp(command, "help"))
  {
    model_usage(model, stderr);
    exit(0);
  }

  // Have we been asked to run a benchmark?
  if (!strcmp(command, "benchmark"))
  {
    if (input == NULL)
    {
      fprintf(stderr, "polymec: No benchmark problem given! Usage:\n");
      fprintf(stderr, "polymec %s benchmark [problem]\n", model_name);
      exit(-1);
    }
    model_run_benchmark(model, input, opts);
    exit(0);
  }

  // We are asked to run a simulation.
  ASSERT(!strcmp(command, "run"));
  if (input == NULL)
  {
    fprintf(stderr, "polymec: No input file given! Usage:\n");
    fprintf(stderr, "polymec %s run [input file]\n", model_name);
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
  options_free(opts);
#endif

  // That's it.
  return 0;
}

#ifdef __cplusplus
}
#endif
