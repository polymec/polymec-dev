// Copyright (c) 2012-2013, Jeffrey N. Johnson
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this 
// list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice, 
// this list of conditions and the following disclaimer in the documentation 
// and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "core/polymec.h"
#include "core/model.h"
#include "core/unordered_map.h"
#include "core/options.h"
#include "core/array.h"
#include "core/array_utils.h"

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

  int save_every; // Save frequency.
  int plot_every; // Plot frequency.

  // Observations.
  string_ptr_unordered_map_t* point_obs;
  string_ptr_unordered_map_t* integ_obs;
  string_ptr_unordered_map_t* min_obs;
  string_ptr_unordered_map_t* max_obs;
  string_array_t* observations;
  real_t* obs_times;
  int num_obs_times, obs_time_index;

  // Data related to a given simulation.
  char* sim_name; // Simulation name.
  real_t time;    // Current simulation time.
  int step;       // Current simulation step.
  real_t max_dt;  // Maximum time step.

  // Interpreter for parsing input files.
  interpreter_t* interpreter;
};

model_t* model_new(const char* name, void* context, model_vtable vtable, options_t* options)
{
  ASSERT(options != NULL);

  model_t* model = malloc(sizeof(model_t));
  model->vtable = vtable;
  model->context = context;
  model->name = string_dup(name);
  model->benchmarks = model_benchmark_map_new();
  model->sim_name = NULL;
  model->save_every = -1;
  model->plot_every = -1;
  model->max_dt = FLT_MAX;
  model->interpreter = NULL;

  // Initialize observation arrays.
  model->point_obs = string_ptr_unordered_map_new();
  model->integ_obs = string_ptr_unordered_map_new();
  model->min_obs = string_ptr_unordered_map_new();
  model->max_obs = string_ptr_unordered_map_new();
  model->observations = string_array_new();
  model->obs_times = NULL;
  model->obs_time_index = 0;

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

  if (model->interpreter != NULL)
    interpreter_free(model->interpreter);

  // Clear observations.
  free(model->obs_times);
  string_array_free(model->observations);
  string_ptr_unordered_map_free(model->min_obs);
  string_ptr_unordered_map_free(model->max_obs);
  string_ptr_unordered_map_free(model->integ_obs);
  string_ptr_unordered_map_free(model->point_obs);

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
  metadata->description = string_dup(description);
  model_benchmark_map_insert_with_kv_dtor(model->benchmarks, string_dup(benchmark), metadata, free_benchmark_kv);
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
      real_t expected_rate = atof(exp_conv_rate);
      real_t actual_rate = atof(conv_rate);
      real_t actual_rate_sigma = atof(sigma);
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
        real_t expected_norm = atof(exp_err_norm);
        real_t actual_norm = atof(err_norm);
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

typedef struct 
{
  real_t (*func)(void*, point_t*, real_t);
  point_t point;
} point_obs_t;

static void point_obs_dtor(char* key, void* val)
{
  free(key);
  point_obs_t* data = val;
  free(data);
}

void model_define_point_observation(model_t* model, 
                                    const char* name, 
                                    real_t (*compute_point_observation)(void* context, 
                                                                        point_t* x,
                                                                        real_t t),
                                    point_t* point)

{
  point_obs_t* data = malloc(sizeof(point_obs_t));
  data->func = compute_point_observation;
  data->point = *point;
  string_ptr_unordered_map_insert_with_kv_dtor(model->point_obs, 
                                               string_dup(name), 
                                               data, point_obs_dtor);
}

static void key_dtor(char* key)
{
  free(key);
}

void model_define_integrated_observation(model_t* model, 
                                         const char* name, 
                                         real_t (*compute_integrated_observation)(void* context,
                                                                                  real_t t))
{
  string_ptr_unordered_map_insert_with_k_dtor(model->integ_obs, 
                                              string_dup(name), 
                                              compute_integrated_observation, 
                                              key_dtor);
}

void model_define_max_observation(model_t* model, 
                                  const char* name, 
                                  real_t (*compute_max_observation)(void* context, 
                                                                    real_t t))
{
  string_ptr_unordered_map_insert_with_k_dtor(model->max_obs, 
                                              string_dup(name), 
                                              compute_max_observation,
                                              key_dtor);
}

void model_define_min_observation(model_t* model, 
                                  const char* name, 
                                  real_t (*compute_min_observation)(void* context, 
                                                                    real_t t))
{
  string_ptr_unordered_map_insert_with_k_dtor(model->min_obs, 
                                              string_dup(name), 
                                              compute_min_observation,
                                              key_dtor);
}

void model_observe(model_t* model, const char* observation)
{
  // Make sure this is an observation in the model.
  if (!string_ptr_unordered_map_contains(model->point_obs, (char*)observation) && 
      !string_ptr_unordered_map_contains(model->integ_obs, (char*)observation) && 
      !string_ptr_unordered_map_contains(model->min_obs, (char*)observation) && 
      !string_ptr_unordered_map_contains(model->max_obs, (char*)observation))
  {
    polymec_error("Observation '%s' is not defined for this model.", observation);
  }
  string_array_append_with_dtor(model->observations, string_dup(observation), key_dtor);
}

void model_set_observation_times(model_t* model, real_t* times, int num_times)
{
  ASSERT(num_times > 0);

  if (model->obs_times != NULL)
    free(model->obs_times);
  model->num_obs_times = num_times;
  model->obs_times = malloc(sizeof(real_t) * num_times);
  memcpy(model->obs_times, times, sizeof(real_t) * num_times);

  // Sort the observation times.
  real_qsort(model->obs_times, model->num_obs_times);
}

static void model_do_periodic_work(model_t* model)
{
  // Do plots and saves.
  if ((model->plot_every > 0) && (model->step % model->plot_every) == 0)
    model_plot(model);
  if ((model->save_every > 0) && (model->step % model->save_every) == 0)
    model_save(model);

  // Now record any observations we need to.
  if (model->time >= model->obs_times[model->obs_time_index])
  {
    model_record_observations(model);
    ++model->obs_time_index;
  }
}

// Initialize the model at the given time.
void model_init(model_t* model, real_t t)
{
  log_detail("%s: Initializing at time %g.", model->name, t);
  model->vtable.init(model->context, t);
  model->step = 0;
  model->time = t;
  model->obs_time_index = 0;

  model_do_periodic_work(model);
}

// Returns the largest permissible time step that can be taken by the model
// starting at time t.
real_t model_max_dt(model_t* model, char* reason)
{
  real_t dt = FLT_MAX;
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

void model_advance(model_t* model, real_t dt)
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
  if (model->sim_name == NULL)
    polymec_error("No simulation name was set with model_set_sim_name.");
  char prefix[strlen(model->sim_name) + 16];
  snprintf(prefix, strlen(model->sim_name) + 16, "%s-%d", model->sim_name, step);
  log_detail("%s: Loading save file from directory %s...", model->name, model->sim_name);
  model->vtable.load(model->context, prefix, model->sim_name, &model->time, step);
  model->step = step;
}

void model_save(model_t* model)
{
  if (model->sim_name == NULL)
    polymec_error("No simulation name was set with model_set_sim_name.");
  char prefix[strlen(model->sim_name) + 16];
  snprintf(prefix, strlen(model->sim_name) + 16, "%s-%d", model->sim_name, model->step);
  log_detail("%s: Writing save file to directory %s...", model->name, model->sim_name);
  model->vtable.save(model->context, prefix, model->sim_name, model->time, model->step);
}

void model_plot(model_t* model)
{
  if (model->sim_name == NULL)
    polymec_error("No simulation name was set with model_set_sim_name.");
  char prefix[strlen(model->sim_name) + 16];
  snprintf(prefix, strlen(model->sim_name) + 16, "%s-%d", model->sim_name, model->step);
  log_detail("%s: Writing plot to directory %s...", model->name, model->sim_name);
  model->vtable.plot(model->context, prefix, model->sim_name, model->time, model->step);
}

void model_record_observations(model_t* model)
{
  if (model->sim_name == NULL)
    polymec_error("No simulation name was set with model_set_sim_name.");

  if (model->observations->size == 0)
    return;

  // Open up the observation file.
  char obs_fn[strlen(model->sim_name) + 5];
  snprintf(obs_fn, strlen(model->sim_name) + 4, "%s.obs", model->sim_name);
  FILE* obs = fopen(obs_fn, "r");
  bool first_time = (obs == NULL);
  if (!first_time)
    fclose(obs);
  obs = fopen(obs_fn, "w+");

  // If we're writing the first entry, write a header
  if (first_time)
  {
    fprintf(obs, "# Observations for %s\n", model->sim_name);
    fprintf(obs, "# time ");
    for (int i = 0; i < model->observations->size; ++i)
      fprintf(obs, "%s ", model->observations->data[i]);
    fprintf(obs, "\n");
  }

  // Write the current simulation time.
  fprintf(obs, "%g ", model->time);

  // Go through all the desired observations.
  for (int i = 0; i < model->observations->size; ++i)
  {
    char* obs_name = model->observations->data[i];

    // Is this one a point observation?
    point_obs_t** point_obs_data = (point_obs_t**)string_ptr_unordered_map_get(model->point_obs, obs_name);
    double value;
    if (point_obs_data != NULL)
      value = (*point_obs_data)->func(model->context, &((*point_obs_data)->point), model->time);
    else
    {
      // Is it an integrated observation?
      typedef real_t (*integ_observation_func)(void*, real_t);
      integ_observation_func* integ_func = (integ_observation_func*)string_ptr_unordered_map_get(model->integ_obs, obs_name);
      if (integ_func != NULL)
        value = (*integ_func)(model->context, model->time);
      else
      {
        // Is it an minimum value observation?
        typedef real_t (*min_observation_func)(void*, real_t);
        min_observation_func* min_func = (min_observation_func*)string_ptr_unordered_map_get(model->min_obs, obs_name);
        if (min_func != NULL)
          value = (*min_func)(model->context, model->time);
        else
        {
          // It's an minimum value observation.
          typedef real_t (*max_observation_func)(void*, real_t);
          max_observation_func* max_func = (max_observation_func*)string_ptr_unordered_map_get(model->min_obs, obs_name);
          ASSERT(max_func != NULL);
          value = (*max_func)(model->context, model->time);
        }
      }
    }

    // Write the observation to the file.
    fprintf(obs, "%g ", value);
  }

  // Close the file.
  fclose(obs);
}

void model_compute_error_norms(model_t* model, st_func_t* solution, real_t* error_norms)
{ 
  if (model->vtable.compute_error_norms != NULL)
    model->vtable.compute_error_norms(model->context, solution, model->time, error_norms);
  else
    polymec_error("%s: This model is not equipped to compute error norms.", model->name);
}

void model_run(model_t* model, real_t t1, real_t t2, int max_steps)
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
      real_t dt = model_max_dt(model, reason);
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
  model->sim_name = string_dup(sim_name);
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
  real_t t1 = 0.0, t2 = 1.0;
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

static void minimal_driver_usage(const char* model_name, FILE* stream)
{
  fprintf(stream, "%s: usage:\n", model_name);
  fprintf(stream, "%s input [option1=val [option2=val2 [...]]]\n\n", model_name);
  exit(-1);
}

int model_minimal_main(const char* model_name, model_ctor constructor, int argc, char* argv[])
{
  // Start everything up.
  polymec_init(argc, argv);

  // Parse options on the command line.
  options_t* opts = options_parse(argc, argv);
  if (opts == NULL)
    minimal_driver_usage(model_name, stderr);

  // Here, the command serves as the input.
  char* input = options_command(opts);
  if ((input == NULL) || !strcmp(input, "help"))
  {
    minimal_driver_usage(model_name, stderr);
  }

  // Attempt to construct the model.
  model_t* model = (*constructor)(opts);
  ASSERT(model != NULL);

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
  real_t t1 = 0.0, t2 = 1.0;
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

void model_report_conv_rate(options_t* options, real_t conv_rate, real_t sigma)
{
  ASSERT(options != NULL);
  char crstr[1024], sigmastr[1024];
  snprintf(crstr, 1024, "%g", conv_rate);
  snprintf(sigmastr, 1024, "%g", sigma);
  options_set(options, "conv_rate", crstr);
  options_set(options, "conv_rate_sigma", sigmastr);
}

void model_report_error_norm(options_t* options, real_t error_norm)
{
  ASSERT(options != NULL);
  char nstr[1024];
  snprintf(nstr, 1024, "%g", error_norm);
  options_set(options, "error_norm", nstr);
}

