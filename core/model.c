#include <stdlib.h>
#include <string.h>
#include <float.h>
#include "core/model.h"
#include "core/unordered_map.h"
#include "core/options.h"

#ifdef __cplusplus
extern "C" {
#endif


// Benchmark metadatum.
typedef struct
{
  model_benchmark_function_t function;
  char* description;
} model_benchmark_t;

// Destructor for benchmark key/value pairs.
static void free_benchmark_kv(char* key, model_benchmark_t* value)
{
  free(key);
  free(value->description);
  free(value);
}

// Mapping from benchmark names to metadata.
DEFINE_UNORDERED_MAP(model_benchmark_map, char*, model_benchmark_t*, string_hash, string_equals)

struct model_t 
{
  // Model metadata.
  void* context;
  char* name;
  model_vtable vtable;
  model_benchmark_map_t* benchmarks;

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

  // Interpreter for parsing input files.
  interpreter_t* interpreter;
};

model_t* model_new(const char* name, void* context, model_vtable vtable, options_t* options)
{
  ASSERT(options != NULL);

  model_t* model = malloc(sizeof(model_t));
  model->vtable = vtable;
  model->context = context;
  model->name = strdup(name);
  model->benchmarks = model_benchmark_map_new();
  model->sim_name = NULL;
  model->saver = NULL;
  model->save_every = -1;
  model->plotter = NULL;
  model->plot_every = -1;
  model->max_dt = FLT_MAX;
  model->interpreter = NULL;

  // Some generic options.
  char* logging = options_value(options, "logging");
  if (logging != NULL)
  {
    if (!strcasecmp(logging, "debug"))
      set_log_level(LOG_DEBUG);
    else if (!strcasecmp(logging, "detail"))
      set_log_level(LOG_DETAIL);
    else if (!strcasecmp(logging, "info"))
      set_log_level(LOG_INFO);
    else if (!strcasecmp(logging, "urgent"))
      set_log_level(LOG_URGENT);
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
  model_benchmark_map_free(model->benchmarks);

  if (model->sim_name != NULL)
    free(model->sim_name);
  if (model->saver != NULL)
    io_free(model->saver);
  if (model->plotter != NULL)
    io_free(model->plotter);

  if (model->interpreter != NULL)
    interpreter_free(model->interpreter);

  free(model);
}

char* model_name(model_t* model)
{
  return model->name;
}

interpreter_t* model_interpreter(model_t* model)
{
  if (model->interpreter == NULL)
    model_enable_interpreter(model, NULL);
  ASSERT(model->interpreter != NULL);
  return model->interpreter;
}

void model_enable_interpreter(model_t* model, interpreter_validation_t* valid_inputs)
{
  if (model->interpreter != NULL)
    interpreter_free(model->interpreter);
  model->interpreter = interpreter_new(valid_inputs);
}

void model_register_benchmark(model_t* model, const char* benchmark, model_benchmark_function_t function, const char* description)
{
  ASSERT(benchmark != NULL);
  ASSERT(function != NULL);
  model_benchmark_t* metadata = malloc(sizeof(model_benchmark_t));
  metadata->function = function;
  metadata->description = strdup(description);
  model_benchmark_map_insert_with_kv_dtor(model->benchmarks, strdup(benchmark), metadata, free_benchmark_kv);
}

void model_run_all_benchmarks(model_t* model, options_t* options)
{
  int pos = 0;
  char *benchmark;
  model_benchmark_t* metadata;
  while (model_benchmark_map_next(model->benchmarks, &pos, &benchmark, &metadata))
    (*metadata->function)(options);
}

void model_run_benchmark(model_t* model, const char* benchmark, options_t* options)
{
  // Try to retrieve this benchmark.
  model_benchmark_t** metadata = model_benchmark_map_get(model->benchmarks, (char*)benchmark);
  if (metadata != NULL)
  {
    // By default (unless overridden), benchmarks communicate only with 
    // "urgent" log messages--others are suppressed.
    if (options_value(options, "logging") == NULL)
      set_log_level(LOG_URGENT);

    log_info("%s: Running benchmark '%s'.", model->name, benchmark);
    options_set(options, "sim_name", benchmark);
    (*(*metadata)->function)(options);
    log_info("%s: Finished running benchmark '%s'.", model->name, benchmark);

    char* conv_rate = options_value(options, "conv_rate");
    char* sigma = options_value(options, "conv_rate_sigma");
    char* exp_conv_rate = options_value(options, "expected_conv_rate");
    if ((exp_conv_rate != NULL) && (conv_rate != NULL))
    {
      double expected_rate = atof(exp_conv_rate);
      double actual_rate = atof(conv_rate);
      double actual_rate_sigma = atof(sigma);
      log_urgent("%s: Expected convergence rate: %g", model->name, expected_rate);
      if (actual_rate_sigma != 0.0)
        log_urgent("%s: Measured convergence rate: %g +/- %g", model->name, actual_rate, actual_rate_sigma);
      else
        log_urgent("%s: Measured convergence rate: %g", model->name, actual_rate);
      if (actual_rate >= expected_rate)
        log_urgent("%s: Benchmark '%s' convergence test PASSED\n", model->name, benchmark);
      else
        log_urgent("%s: Benchmark '%s' convergence test FAILED\n", model->name, benchmark);
    }
    else
    {
      char* exp_err_norm = options_value(options, "expected_error_norm");
      char* err_norm = options_value(options, "error_norm");
      if ((exp_err_norm != NULL) && (err_norm != NULL))
      {
        double expected_norm = atof(exp_err_norm);
        double actual_norm = atof(err_norm);
        log_urgent("%s: Expected error norm: %g", model->name, expected_norm);
        log_urgent("%s: Measured error norm: %g", model->name, actual_norm);
        if (actual_norm <= expected_norm)
          log_urgent("%s: Benchmark '%s' error norm test PASSED\n", model->name, benchmark);
        else
          log_urgent("%s: Benchmark '%s' error norm test FAILED\n", model->name, benchmark);
      }
    }
  }
  else
  {
    polymec_error("%s: Benchmark not found: '%s'.", model->name, benchmark);
  }
}

void model_read_input_string(model_t* model, const char* input, options_t* options)
{
  interpreter_t* interp = model_interpreter(model);
  interpreter_parse_string(interp, (char*)input);

  // Load the inputs into the model.
  model->vtable.read_input(model->context, interp, options);
}

void model_read_input_file(model_t* model, const char* file, options_t* options)
{
  interpreter_t* interp = model_interpreter(model);
  interpreter_parse_file(interp, (char*)file);

  bool no_opts = false;
  if (options == NULL)
  {
    no_opts = true;
    options = options_new();
  }

  // Load the inputs into the model.
  model->vtable.read_input(model->context, interp, options);

  if (no_opts)
    options = NULL; 
}

static void model_do_periodic_work(model_t* model)
{
  if ((model->plot_every > 0) && (model->step % model->plot_every) == 0)
    model_plot(model);
  if ((model->save_every > 0) && (model->step % model->save_every) == 0)
    model_save(model);
}

// Initialize the model at the given time.
void model_init(model_t* model, double t)
{
  log_detail("%s: Initializing at time %g.", model->name, t);
  model->vtable.init(model->context, t);
  model->step = 0;
  model->time = t;

  model_do_periodic_work(model);
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
  log_info("%s: Step %d (t = %g, dt = %g)", model->name, model->step, model->time, dt);
  model->vtable.advance(model->context, model->time, dt);
  model->time += dt;
  model->step += 1;

  model_do_periodic_work(model);
}

void model_finalize(model_t* model)
{
  log_detail("%s: Finalizing model at t = %g", model->name, model->time);
  if (model->vtable.finalize != NULL)
    model->vtable.finalize(model->context, model->step, model->time);
}

void model_load(model_t* model, int step)
{
  ASSERT(step >= 0);
  if (model->saver == NULL)
    polymec_error("No saver/loader was set with model_set_saver.");
  if (model->sim_name == NULL)
    polymec_error("No simulation name was set with model_set_sim_name.");
  char prefix[strlen(model->sim_name) + 16];
  snprintf(prefix, strlen(model->sim_name) + 16, "%s-%d", model->sim_name, step);
  log_detail("%s: Loading save file from directory %s...", model->name, model->sim_name);
  io_open(model->saver, prefix, model->sim_name, IO_READ);
  model->vtable.load(model->context, model->saver, &model->time, step);
  model->step = step;
  io_close(model->saver);
}

void model_save(model_t* model)
{
  if (model->saver == NULL)
    polymec_error("No saver/loader was set with model_set_saver.");
  if (model->sim_name == NULL)
    polymec_error("No simulation name was set with model_set_sim_name.");
  char prefix[strlen(model->sim_name) + 16];
  snprintf(prefix, strlen(model->sim_name) + 16, "%s-%d", model->sim_name, model->step);
  log_detail("%s: Writing save file to directory %s...", model->name, model->sim_name);
  io_open(model->saver, prefix, model->sim_name, IO_WRITE);
  model->vtable.save(model->context, model->saver, model->time, model->step);
  io_close(model->saver);
}

void model_plot(model_t* model)
{
  if (model->plotter == NULL)
    polymec_error("No plotter was set with model_set_plotter.");
  if (model->sim_name == NULL)
    polymec_error("No simulation name was set with model_set_sim_name.");
  char prefix[strlen(model->sim_name) + 16];
  snprintf(prefix, strlen(model->sim_name) + 16, "%s-%d", model->sim_name, model->step);
  log_detail("%s: Writing plot to directory %s...", model->name, model->sim_name);
  io_open(model->plotter, prefix, model->sim_name, IO_WRITE);
  model->vtable.plot(model->context, model->plotter, model->time, model->step);
  io_close(model->plotter);
}

void model_compute_error_norms(model_t* model, st_func_t* solution, double* error_norms)
{ 
  if (model->vtable.compute_error_norms != NULL)
    model->vtable.compute_error_norms(model->context, solution, model->time, error_norms);
  else
    polymec_error("%s: This model is not equipped to compute error norms.", model->name);
}

void model_run(model_t* model, double t1, double t2, int max_steps)
{
  ASSERT(t2 >= t1);
  if (t2 > t1)
  {
    if (max_steps == INT_MAX)
      log_detail("%s: Running from time %g to %g.", model->name, t1, t2);
    else
      log_detail("%s: Running from time %g to %g (or for %d steps).", model->name, t1, t2, max_steps);
  }
  else
    log_detail("%s: Running simulation at time %g.", model->name, t1);
  model_init(model, t1);

  if (t1 == t2)
  {
    model_do_periodic_work(model);
  }
  else
  {
    while ((model->time < t2) && (model->step < max_steps))
    {
      char reason[POLYMEC_MODEL_MAXDT_REASON_SIZE];
      double dt = model_max_dt(model, reason);
      if (dt > t2 - model->time)
      {
        dt = t2 - model->time;
        snprintf(reason, POLYMEC_MODEL_MAXDT_REASON_SIZE, "End of simulation");
      }
      log_detail("%s: Selected time step dt = %g\n (Reason: %s).", model->name, dt, reason);
      model_advance(model, dt);
    }
  }

  // Do any finalization.
  model_finalize(model);

  log_detail("%s: Run concluded at time %g.", model->name, t2);
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
  if (model->saver != NULL)
    io_free(model->saver);
  model->saver = saver;
}

void model_set_plotter(model_t* model, io_interface_t* plotter)
{
  ASSERT(plotter != NULL);
  ASSERT(plotter != model->saver);
  if (model->plotter != NULL)
    io_free(model->plotter);
  model->plotter = plotter;
}

int model_main(const char* model_name, model_ctor constructor, int argc, char* argv[])
{
  // Start everything up.
  polymec_init(argc, argv);

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
  static const char* valid_commands[] = {"run", "benchmark", "list-benchmarks", "help", NULL};
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
    char *benchmark;
    model_benchmark_t* metadata;
    while (model_benchmark_map_next(model->benchmarks, &pos, &benchmark, &metadata))
      fprintf(stderr, "  %s (%s)\n", benchmark, metadata->description);
    fprintf(stderr, "\n");
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

  // By default, the simulation is named after the input file (minus its suffix).
  model_set_sim_name(model, input);

  // Read the contents of the input file into the model's interpreter.
  model_read_input_file(model, input, opts);

  // Default time endpoints, max number of steps.
  double t1 = 0.0, t2 = 1.0;
  int max_steps = INT_MAX;

  // Overwrite these defaults with interpreted values.
  interpreter_t* interp = model_interpreter(model);
  if (interpreter_contains(interp, "t1", INTERPRETER_NUMBER))
    t1 = interpreter_get_number(interp, "t1");
  if (interpreter_contains(interp, "t2", INTERPRETER_NUMBER))
    t2 = interpreter_get_number(interp, "t2");
  if (interpreter_contains(interp, "max_steps", INTERPRETER_NUMBER))
    max_steps = (int)interpreter_get_number(interp, "max_steps");

  // If these are given as options, they are overridden by the command line.
  {
    char* opt = options_value(opts, "t1");
    if (opt != NULL)
      t1 = atof(opt);
    opt = options_value(opts, "t2");
    if (opt != NULL)
      t2 = atof(opt);
    opt = options_value(opts, "max_steps");
    if (opt != NULL)
      max_steps = atoi(opt);
  }

  // Run the model.
  model_run(model, t1, t2, max_steps);

  // Clean up.
  model_free(model);

  return 0;
}

void model_report_conv_rate(options_t* options, double conv_rate, double sigma)
{
  ASSERT(options != NULL);
  char crstr[1024], sigmastr[1024];
  snprintf(crstr, 1024, "%g", conv_rate);
  snprintf(sigmastr, 1024, "%g", sigma);
  options_set(options, "conv_rate", crstr);
  options_set(options, "conv_rate_sigma", sigmastr);
}

void model_report_error_norm(options_t* options, double error_norm)
{
  ASSERT(options != NULL);
  char nstr[1024];
  snprintf(nstr, 1024, "%g", error_norm);
  options_set(options, "error_norm", nstr);
}

#ifdef __cplusplus
}
#endif

