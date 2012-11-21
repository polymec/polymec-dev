#include <stdlib.h>
#include <string.h>
#include <float.h>
#include "core/model.h"
#include "core/unordered_map.h"
#include "core/options.h"

#ifdef __cplusplus
extern "C" {
#endif

struct model_t 
{
  // Model metadata.
  void* context;
  char* name;
  model_vtable vtable;
  str_str_unordered_map_t* benchmarks;

  // Functions that are called periodically. 
  io_interface_t* saver; // I/O interface for saving.
  int save_every; // Save frequency.
  io_interface_t* plotter; // I/O interface for plotting.
  int plot_every; // Plot frequency.

  // Data related to a given simulation.
  char* sim_name; // Simulation name.
  double time;    // Current simulation time.
  int step;       // Current simulation step.
  double max_dt;  // Maximum time step.
};

model_t* model_new(const char* name, void* context, model_vtable vtable, options_t* options)
{
  ASSERT(options != NULL);

  model_t* model = malloc(sizeof(model_t));
  model->vtable = vtable;
  model->context = context;
  model->name = strdup(name);
  model->benchmarks = str_str_unordered_map_new();
  model->sim_name = NULL;
  model->saver = NULL;
  model->save_every = -1;
  model->plotter = NULL;
  model->plot_every = -1;
  model->max_dt = FLT_MAX;

  // Some generic options.
  char* logging = options_value(options, "logging");
  if (logging != NULL)
  {
    if (!strcasecmp(logging, "debug"))
      set_log_level(LOG_DEBUG);
    else if (!strcasecmp(logging, "info"))
      set_log_level(LOG_INFO);
    else if (!strcasecmp(logging, "warning"))
      set_log_level(LOG_WARNING);
    else if (!strcasecmp(logging, "off"))
      set_log_level(LOG_NONE);
  }
  char* plot_every = options_value(options, "plot_every");
  if (plot_every != NULL)
    model->plot_every = atoi(plot_every);
  char* save_every = options_value(options, "save_every");
  if (save_every != NULL)
    model->save_every = atoi(save_every);
  char* max_dt = options_value(options, "max_dt");
  if (max_dt != NULL)
    model->max_dt = atof(max_dt);
  char* sim_name = options_value(options, "sim_name");
  if (sim_name != NULL)
    model_set_sim_name(model, sim_name);
  return model;
}

void model_free(model_t* model)
{
  if ((model->context != NULL) && (model->vtable.dtor != NULL))
    model->vtable.dtor(model->context);
  free(model->name);

  // Clear benchmarks.
  str_str_unordered_map_free(model->benchmarks);

  if (model->sim_name != NULL)
    free(model->sim_name);
  if (model->saver != NULL)
    io_free(model->saver);
  if (model->plotter != NULL)
    io_free(model->plotter);
  free(model);
}

char* model_name(model_t* model)
{
  return model->name;
}

static void destroy_key_and_value(char* key, char* value)
{
  free(key);
  free(value);
}

void model_register_benchmark(model_t* model, const char* benchmark, const char* description)
{
  str_str_unordered_map_insert_with_dtor(model->benchmarks, strdup(benchmark), strdup(description), destroy_key_and_value);
}

void model_run_all_benchmarks(model_t* model, options_t* options)
{
  int pos = 0;
  char *benchmark, *descr;
  while (str_str_unordered_map_next(model->benchmarks, &pos, &benchmark, &descr))
    model_run_benchmark(model, benchmark, options);
}

void model_run_benchmark(model_t* model, const char* benchmark, options_t* options)
{
  if (model->vtable.run_benchmark != NULL)
  {
    log_info("%s: Running benchmark '%s'.", model->name, benchmark);
    options_set(options, "sim_name", benchmark);
    model->vtable.run_benchmark(benchmark, options);
    log_info("%s: Finished running benchmark '%s'.", model->name, benchmark);
  }
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
  log_info("%s: Initializing at time %g.", model->name, t);
  model->vtable.init(model->context, t);
  model->step = 0;
  model->time = t;
}

// Returns the largest permissible time step that can be taken by the model
// starting at time t.
double model_max_dt(model_t* model, char* reason)
{
  double dt = FLT_MAX;
  strcpy(reason, "No time step constraints.");
  if (model->max_dt < FLT_MAX)
  {
    dt = model->max_dt;
    strcpy(reason, "Max dt set in options.");
  }
  if (model->vtable.max_dt != NULL)
    return model->vtable.max_dt(model->context, model->time, reason);
  return dt;
}

void model_advance(model_t* model, double dt)
{
  log_info("%s: Advancing from time %g to %g.", model->name, model->time, model->time+dt);
  model->vtable.advance(model->context, model->time, dt);
  model->time += dt;
  model->step += 1;

  // Call periodic work.
  if ((model->plot_every > 0) && (model->plot_every % model->step) == 0)
    model_plot(model);
  if ((model->save_every > 0) && (model->save_every % model->step) == 0)
    model_save(model);
}

void model_load(model_t* model, int step)
{
  ASSERT(step >= 0);
  if (model->saver == NULL)
    arbi_error("No saver/loader was set with model_set_saver.");
  if (model->sim_name == NULL)
    arbi_error("No simulation name was set with model_set_sim_name.");
  char prefix[strlen(model->sim_name) + 16];
  snprintf(prefix, strlen(model->sim_name) + 16, "%s-%d", model->sim_name, step);
  io_open(model->saver, prefix, model->sim_name, IO_READ);
  model->vtable.load(model->context, model->saver, &model->time, step);
  model->step = step;
  io_close(model->saver);
}

void model_save(model_t* model)
{
  if (model->saver == NULL)
    arbi_error("No saver/loader was set with model_set_saver.");
  if (model->sim_name == NULL)
    arbi_error("No simulation name was set with model_set_sim_name.");
  char prefix[strlen(model->sim_name) + 16];
  snprintf(prefix, strlen(model->sim_name) + 16, "%s-%d", model->sim_name, model->step);
  io_open(model->saver, prefix, model->sim_name, IO_WRITE);
  model->vtable.save(model->context, model->saver, model->time, model->step);
  io_close(model->saver);
}

void model_plot(model_t* model)
{
  if (model->plotter == NULL)
    arbi_error("No plotter was set with model_set_plotter.");
  if (model->sim_name == NULL)
    arbi_error("No simulation name was set with model_set_sim_name.");
  char prefix[strlen(model->sim_name) + 16];
  snprintf(prefix, strlen(model->sim_name) + 16, "%s-%d", model->sim_name, model->step);
  io_open(model->plotter, prefix, model->sim_name, IO_WRITE);
  model->vtable.plot(model->context, model->plotter, model->time, model->step);
  io_close(model->plotter);
}

void model_run(model_t* model, double t1, double t2)
{
  log_info("%s: Running from time %g to %g.", model->name, t1, t2);
  model_init(model, t1);
  while (model->time < t2)
  {
    char reason[ARBI_MODEL_MAXDT_REASON_SIZE];
    double dt = model_max_dt(model, reason);
    if (dt > t2 - model->time)
    {
      dt = t2 - model->time;
      snprintf(reason, ARBI_MODEL_MAXDT_REASON_SIZE, "End of simulation");
    }
    log_info("%s: Selected time step dt = %g\n (Reason: %s).", model->name, dt, reason);
    model_advance(model, dt);
  }
  log_info("%s: Run concluded at time %g.", model->name, t2);
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

void model_set_sim_name(model_t* model, const char* sim_name)
{
  ASSERT(sim_name != NULL);
  if (model->sim_name != NULL)
    free(model->sim_name);
  model->sim_name = strdup(sim_name);
}

void model_set_saver(model_t* model, io_interface_t* saver)
{
  ASSERT(saver != NULL);
  ASSERT(saver != model->plotter);
  model->saver = saver;
}

void model_set_plotter(model_t* model, io_interface_t* plotter)
{
  ASSERT(plotter != NULL);
  ASSERT(plotter != model->saver);
  model->plotter = plotter;
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
      model_run_all_benchmarks(model, opts);
    else
      model_run_benchmark(model, input, opts);

    model_free(model);
    return 0;
  }

  // Have we been asked to list all available benchmarks?
  if (!strcmp(command, "list-benchmarks"))
  {
    fprintf(stderr, "Benchmarks for %s model:\n", model_name);
    int pos = 0;
    char *benchmark, *descr;
    while (str_str_unordered_map_next(model->benchmarks, &pos, &benchmark, &descr))
      fprintf(stderr, "  %s (%s)\n", benchmark, descr);
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

  // By default, the simulation is named after the input file.
  model_set_sim_name(model, input);

  // Default time endpoints.
  double t1 = 0.0, t2 = 1.0;

  // Run the model.
  model_run(model, t1, t2);

  // Clean up.
  model_free(model);

  return 0;
}

#ifdef __cplusplus
}
#endif

