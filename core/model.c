#include <stdlib.h>
#include <string.h>
#include <float.h>
#include "core/model.h"
#include "core/options.h"
#include "core/simulation.h"

#ifdef __cplusplus
extern "C" {
#endif

struct model_t 
{
  void* context;
  char* name;
  model_vtable vtable;
  char** benchmarks;
  int num_benchmarks;
};

model_t* model_new(const char* name, void* context, model_vtable vtable)
{
  model_t* model = malloc(sizeof(model_t));
  model->vtable = vtable;
  model->context = context;
  model->name = strdup(name);
  model->num_benchmarks = 0;
  model->benchmarks = NULL;
  return model;
}

void model_free(model_t* model)
{
  if ((model->context != NULL) && (model->vtable.dtor != NULL))
    model->vtable.dtor(model->context);
  free(model->name);
  model_register_benchmarks(model, NULL);
  free(model);
}

char* model_name(model_t* model)
{
  return model->name;
}

void model_register_benchmarks(model_t* model, const char** benchmarks)
{
  if (model->benchmarks != NULL)
  {
    for (int i = 0; i < model->num_benchmarks; ++i)
      free(model->benchmarks[i]);
    free(model->benchmarks);
    model->benchmarks = NULL;
  }
  model->num_benchmarks = 0;
  if (benchmarks != NULL)
  {
    char* p = (char*)benchmarks[0];
    while (p != NULL) 
    {
      model->num_benchmarks++;
      p = (char*)benchmarks[model->num_benchmarks];
    }
    model->benchmarks = malloc(sizeof(char*)*model->num_benchmarks);
    for (int i = 0; i < model->num_benchmarks; ++i)
      model->benchmarks[i] = strdup(benchmarks[i]);
  }
}

void model_get_benchmarks(model_t* model, char*** benchmarks, int* num_benchmarks)
{
  *benchmarks = model->benchmarks;
  *num_benchmarks = model->num_benchmarks;
}

void model_run_all_benchmarks(model_t* model)
{
  for (int i = 0; i < model->num_benchmarks; ++i)
    model_run_benchmark(model, (const char*)model->benchmarks[i]);
}

void model_run_benchmark(model_t* model, const char* benchmark)
{
  if (model->vtable.run_benchmark != NULL)
    model->vtable.run_benchmark(benchmark);
  else
  {
    char err[1024];
    snprintf(err, 1024, "No benchmarks are defined for the model '%s'.", model->name);
    arbi_error(err);
  }
}

// Initialize the model at the given time.
void model_init(model_t* model, double t)
{
  model->vtable.init(model->context, t);
}

// Returns the largest permissible time step that can be taken by the model
// starting at time t.
double model_max_dt(model_t* model, double t, char* reason)
{
  if (model->vtable.max_dt != NULL)
    return model->vtable.max_dt(model->context, t, reason);
  else
  {
    strcpy(reason, "No time step constraints.");
    return FLT_MAX;
  }
}

void model_advance(model_t* model, double t, double dt)
{
  model->vtable.advance(model->context, t, dt);
}

void model_load(model_t* model, io_interface_t* io, double* t, int* step)
{
  model->vtable.load(model->context, io, t, step);
}

void model_dump(model_t* model, io_interface_t* io, double t, int step)
{
  model->vtable.dump(model->context, io, t, step);
}

void model_plot(model_t* model, plot_interface_t* plot, double t, int step)
{
  model->vtable.plot(model->context, plot, t, step);
}

void model_run(model_t* model, double t1, double t2)
{
  double t = t1;
  model_init(model, t);
  while (t < t2)
  {
    char reason[ARBI_MODEL_MAXDT_REASON_SIZE];
    double dt = model_max_dt(model, t, reason);
    model_advance(model, t, dt);
    t += dt;
  }
}

void* model_context(model_t* model)
{
  return model->context;
}

static void driver_usage(const char* model_name, FILE* stream)
{
  fprintf(stream, "%s: usage:\n", model_name);
  fprintf(stream, "%s [command] [args]\n\n", model_name);
  fprintf(stream, "Here, [command] is one of the following:\n\n");
  fprintf(stream, "  run [file]           -- Runs a simulation with input from the given file.\n");
  fprintf(stream, "  generate-mesh [file] -- Runs a simulation with input from the given file.\n");
  fprintf(stream, "  benchmark [name]     -- Runs the given benchmark problem ('all' for all).\n");
  fprintf(stream, "  list-benchmarks      -- Lists all benchmark problems.\n");
  fprintf(stream, "  help                 -- Prints information about the given model.\n\n");
  exit(-1);
}

int model_main(const char* model_name, model_ctor constructor, int argc, char* argv[])
{
  // Start everything up.
  arbi_init(argc, argv);

  // Parse options on the command line.
  options_t* opts = options_parse(argc, argv);
  if (opts == NULL)
    driver_usage(model_name, stderr);

  // Extract the command and arguments.
  char* command = options_command(opts);
  char* input = options_input(opts);
  if (!strcmp(command, "help"))
    driver_usage(model_name, stderr);

  // Validate our inputs.
  ASSERT(command != NULL);
  int c = 0;
  static const char* valid_commands[] = {"run", "generate-mesh", "benchmark", "list-benchmarks", "help", NULL};
  while (valid_commands[c] != NULL)
  {
    if (!strcmp(command, valid_commands[c]))
      break;
    ++c;
  }
  if (valid_commands[c] == NULL)
  {
    fprintf(stderr, "%s: invalid command: '%s'\n", model_name, command);
    return -1;
  }

  // Attempt to construct the model.
  model_t* model = (*constructor)(opts);
  ASSERT(model != NULL);

  // Have we been asked to run a benchmark?
  if (!strcmp(command, "benchmark"))
  {
    if (input == NULL)
    {
      fprintf(stderr, "%s: No benchmark problem given! Usage:\n", model_name);
      fprintf(stderr, "%s benchmark [problem]\n", model_name);
      return -1;
    }

    // Have we been asked to run all benchmarks?
    if (!strcmp(input, "all"))
      model_run_all_benchmarks(model);
    else
      model_run_benchmark(model, input);
    return 0;
  }

  // Have we been asked to list all available benchmarks?
  if (!strcmp(command, "list-benchmarks"))
  {
    fprintf(stderr, "Benchmarks for %s model:\n", model_name);
    for (int i = 0; i < model->num_benchmarks; ++i)
      fprintf(stderr, "  %s\n", (const char*)model->benchmarks[i]);
    fprintf(stderr, "\n");
    return 0;
  }

  // Have we been asked to generate a mesh for a problem?
  if (!strcmp(command, "generate-mesh"))
  {
    if (input == NULL)
    {
      fprintf(stderr, "%s: No input file given for mesh generation! Usage:\n", model_name);
      fprintf(stderr, "%s generate-mesh [input file]\n", model_name);
      return -1;
    }

    // Check to see whether the given file exists.
    FILE* fp = fopen(input, "r");
    if (fp == NULL)
    {
      fprintf(stderr, "%s: Input file not found: %s\n", model_name, input);
      return -1;
    }
    fclose(fp);

    // FIXME: Mesh generation goes here!
    return 0;
  }

  // We are asked to run a simulation.
  ASSERT(!strcmp(command, "run"));
  if (input == NULL)
  {
    fprintf(stderr, "%s: No input file given! Usage:\n", model_name);
    fprintf(stderr, "%s run [input file]\n", model_name);
    return -1;
  }

  // Check to see whether the given file exists.
  FILE* fp = fopen(input, "r");
  if (fp == NULL)
  {
    fprintf(stderr, "%s: Input file not found: %s\n", model_name, input);
    return -1;
  }
  fclose(fp);

  // Create a simulation in which to execute the model.
  simulation_t* sim = simulation_new(model, input, opts);

  // Initialize the simulation.
  simulation_init(sim);

  // Run the simulation.
  simulation_run(sim);

  // Clean up.
  simulation_free(sim);
  model_free(model);

  return 0;
}

#ifdef __cplusplus
}
#endif

