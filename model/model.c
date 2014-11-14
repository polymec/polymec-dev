// Copyright (c) 2012-2014, Jeffrey N. Johnson
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

#include <stdarg.h>
#include "core/polymec.h"
#include "core/unordered_map.h"
#include "core/array.h"
#include "core/array_utils.h"
#include "core/text_buffer.h"
#include "model/model.h"

// This helper writes the given message out to stderr on rank 0.
static void print_to_rank0(const char* message, ...)
{
  static int rank = -1;
  if (rank == -1)
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0)
  {
    // Extract the variadic arguments and splat them into a string.
    char m[8192];
    va_list argp;
    va_start(argp, message);
    vsnprintf(m, 8192, message, argp);
    va_end(argp);

    // Write it out.
    fprintf(stderr, "%s", m);
  }
}

// Benchmark metadatum.
typedef struct
{
  model_benchmark_function_t function;
  docstring_t* description;
} model_benchmark_t;

// Destructor for benchmark key/value pairs.
static void free_benchmark_kv(char* key, model_benchmark_t* value)
{
  polymec_free(key);
  polymec_free(value);
}

// Mapping from benchmark names to metadata.
DEFINE_UNORDERED_MAP(model_benchmark_map, char*, model_benchmark_t*, string_hash, string_equals)

struct model_t 
{
  // Model metadata.
  void* context;
  char* name;
  model_vtable vtable;
  docstring_t* doc; // Documentation string.
  model_benchmark_map_t* benchmarks;

  int save_every; // Save frequency.
  int plot_every; // Plot frequency.

  int load_step; // -1 if starting a new sim, >= 0 if loading from a file.

  // Observations.
  string_ptr_unordered_map_t* point_obs;
  string_ptr_unordered_map_t* global_obs;
  string_array_t* observations;
  real_t* obs_times, observe_every;
  int num_obs_times;

  // Data related to a given simulation.
  char* sim_name;    // Simulation name.
  real_t time;       // Current simulation time.
  real_t wall_time0; // Wall time at simulation start.
  real_t wall_time;  // Current wall time.
  real_t sim_speed;  // Simulation "speed" (sim time / wall time per step).
  int step;          // Current simulation step number.
  real_t dt;         // Current simulation time step.
  real_t max_dt;     // Maximum time step.

  // Interpreter for parsing input files.
  interpreter_t* interpreter;
};

// This helper parses a comma-delimited list of observation times, returning 
// an array of times.
static real_t* parse_observation_times(char* observation_time_str, int* num_times)
{
  char** time_strings = string_split(observation_time_str, ",", num_times);
  ASSERT(time_strings != NULL);
  ASSERT(*num_times > 0);
  for (int i = 0; i < *num_times; ++i)
  {
    if (!string_is_number(time_strings[i]))
      polymec_error("Invalid observation time at index %d: %s\n", i, time_strings[i]);
  }
  real_t* times = polymec_malloc(sizeof(real_t) * (*num_times));
  for (int i = 0; i < *num_times; ++i)
  {
    times[i] = (real_t)atof(time_strings[i]);
    polymec_free(time_strings[i]);
  }
  polymec_free(time_strings);

  return times;
}

model_t* model_new(const char* name, void* context, model_vtable vtable, docstring_t* doc)
{
  // First, inspect the virtual table.

  // Exactly one of read_input and read_custom_input should be given.
  ASSERT(((vtable.read_input != NULL) && (vtable.read_custom_input == NULL)) ||
         ((vtable.read_input == NULL) && (vtable.read_custom_input != NULL)));

  model_t* model = polymec_malloc(sizeof(model_t));
  model->vtable = vtable;
  model->context = context;
  model->name = string_dup(name);
  model->doc = doc;
  model->benchmarks = model_benchmark_map_new();
  model->sim_name = NULL;
  model->save_every = -1;
  model->plot_every = -1;
  model->load_step = -1;
  model->observe_every = -FLT_MAX;
  model->wall_time = 0.0;
  model->wall_time0 = 0.0;
  model->dt = 0.0;
  model->step = 0;
  model->max_dt = FLT_MAX;
  model->interpreter = NULL;

  // Initialize observation arrays.
  model->point_obs = string_ptr_unordered_map_new();
  model->global_obs = string_ptr_unordered_map_new();
  model->observations = string_array_new();
  model->num_obs_times = 0;
  model->obs_times = NULL;

  return model;
}

void model_free(model_t* model)
{
  if ((model->context != NULL) && (model->vtable.dtor != NULL))
    model->vtable.dtor(model->context);
  polymec_free(model->name);

  // Clear benchmarks.
  model_benchmark_map_free(model->benchmarks);

  if (model->sim_name != NULL)
    polymec_free(model->sim_name);

  if (model->interpreter != NULL)
    interpreter_free(model->interpreter);

  // Clear observations.
  polymec_free(model->obs_times);
  string_array_free(model->observations);
  string_ptr_unordered_map_free(model->global_obs);
  string_ptr_unordered_map_free(model->point_obs);

  polymec_free(model);
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

void model_register_benchmark(model_t* model, const char* benchmark, model_benchmark_function_t function, docstring_t* description)
{
  ASSERT(benchmark != NULL);
  ASSERT(function != NULL);
  model_benchmark_t* metadata = polymec_malloc(sizeof(model_benchmark_t));
  metadata->function = function;
  metadata->description = description;
  model_benchmark_map_insert_with_kv_dtor(model->benchmarks, string_dup(benchmark), metadata, free_benchmark_kv);
}

void model_describe_benchmark(model_t* model, const char* benchmark, FILE* stream)
{
  model_benchmark_t** metadata_p = (model_benchmark_t**)model_benchmark_map_get(model->benchmarks, (char*)benchmark);
  if (metadata_p != NULL)
  {
    if ((*metadata_p)->description != NULL)
    {
      fprintf(stream, "%s benchmark '%s':\n", model->name, benchmark);
      int pos = 0;
      char* line;
      while (docstring_next((*metadata_p)->description, &pos, &line))
        fprintf(stream, "%s\n", line);
    }
    else
    {
      fprintf(stream, "The benchmark '%s' for the %s model has no description.\n",
              benchmark, model->name);
    }
  }
  else
  {
    fprintf(stream, "No benchmark '%s' was registered for the %s model.\n",
            benchmark, model->name);
  }
}

void model_run_all_benchmarks(model_t* model)
{
  int pos = 0;
  char *benchmark;
  model_benchmark_t* metadata;
  while (model_benchmark_map_next(model->benchmarks, &pos, &benchmark, &metadata))
    (*metadata->function)();
}

void model_run_benchmark(model_t* model, const char* benchmark)
{
  // Try to retrieve this benchmark.
  model_benchmark_t** metadata = model_benchmark_map_get(model->benchmarks, (char*)benchmark);
  options_t* options = options_argv();
  if (metadata != NULL)
  {
    // By default (unless overridden), benchmarks communicate only with 
    // "urgent" log messages--others are suppressed.
    const char* logging = options_value(options, "logging");
    if (logging == NULL)
      set_log_level(LOG_URGENT);
    else if (!strcmp(logging, "debug"))
      set_log_level(LOG_DEBUG);
    else if (!strcmp(logging, "detail"))
      set_log_level(LOG_DETAIL);
    else if (!strcmp(logging, "info"))
      set_log_level(LOG_INFO);

    log_info("%s: Running benchmark '%s'.", model->name, benchmark);
    options_set(options, "sim_name", benchmark);
    (*(*metadata)->function)();
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

static void model_read_input(model_t* model, interpreter_t* interp)
{
  options_t* options = options_argv();

  // We always read certain inputs.
  if (interpreter_contains(interp, "load_step", INTERPRETER_NUMBER))
  {
    model->load_step = (int)interpreter_get_number(interp, "load_step");
    if (model->load_step < 0)
      polymec_error("Invalid load_step: %d (must be non-negative).", model->load_step);
  }
  if (interpreter_contains(interp, "logging", INTERPRETER_STRING))
  {
    char* logging = interpreter_get_string(interp, "logging");
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
    else
      polymec_error("Invalid logging: %s\nMust be one of: debug, detail, info, urgent, off", logging);
  }
  if (interpreter_contains(interp, "plot_every", INTERPRETER_NUMBER))
  {
    model->plot_every = (int)interpreter_get_number(interp, "plot_every");
    if (model->plot_every < 1)
      polymec_error("Invalid (non-positive) plot interval: %d\n", model->plot_every);
  }
  if (interpreter_contains(interp, "save_every", INTERPRETER_NUMBER))
  {
    model->save_every = (int)interpreter_get_number(interp, "save_every");
    if (model->save_every < 1)
      polymec_error("Invalid (non-positive) save interval: %d\n", model->save_every);
  }
  if (interpreter_contains(interp, "observe_every", INTERPRETER_NUMBER))
  {
    model->observe_every = (int)interpreter_get_number(interp, "observe_every");
    if (model->observe_every <= 0.0)
      polymec_error("Invalid (non-positive) observation interval: %g\n", model->observe_every);
  }
  if (interpreter_contains(interp, "observation_times", INTERPRETER_SEQUENCE))
  {
    if (model->observe_every > 0.0)
      polymec_error("Only one of observe_every and observation_times may be specified.");
    int num_obs_times;
    real_t* obs_times = interpreter_get_sequence(interp, "observation_times", &num_obs_times);
    model_set_observation_times(model, obs_times, num_obs_times);
    polymec_free(obs_times);
  }

  // If observation names are given, handle them here.
  if (interpreter_contains(interp, "observations", INTERPRETER_STRING_LIST))
  {
    string_array_clear(model->observations);
    int num_obs;
    char** obs_names = interpreter_get_stringlist(interp, "observations", &num_obs);
    for (int i = 0; i < num_obs; ++i)
    {
      model_observe(model, (const char*)obs_names[i]);
      polymec_free(obs_names[i]);
    }
    polymec_free(obs_names); 
  }

  if (interpreter_contains(interp, "max_dt", INTERPRETER_NUMBER))
  {
    model->max_dt = interpreter_get_number(interp, "max_dt");
    if (model->max_dt <= 0.0)
      polymec_error("Invalid value for max_dt: %g", model->max_dt);
  }

  if (interpreter_contains(interp, "sim_name", INTERPRETER_STRING))
    model_set_sim_name(model, interpreter_get_string(interp, "sim_name"));

  // Read the model-specific inputs.
  model->vtable.read_input(model->context, interp, options);
}

void model_read_input_string(model_t* model, const char* input)
{
  if (model->vtable.read_input != NULL)
  {
    interpreter_t* interp = model_interpreter(model);
    interpreter_parse_string(interp, (char*)input);
    model_read_input(model, interp);
  }
  else
  {
    options_t* options = options_argv();
    model->vtable.read_custom_input(model->context, input, options);
  }
}

void model_read_input_file(model_t* model, const char* file)
{
  if (model->vtable.read_input != NULL)
  {
    interpreter_t* interp = model_interpreter(model);
    interpreter_parse_file(interp, (char*)file);
    model_read_input(model, interp);
  }
  else
  {
    log_detail("%s: Reading custom input from '%s'...", model->name, file);
    text_buffer_t* text = text_buffer_from_file(file);
    char* input = text_buffer_to_string(text);
    text_buffer_free(text);
    options_t* options = options_argv();
    model->vtable.read_custom_input(model->context, (const char*)input, options);
    polymec_free(input);
  }
}

typedef struct 
{
  real_t (*func)(void*, point_t*, real_t);
  point_t point;
} point_obs_t;

static void point_obs_dtor(char* key, void* val)
{
  polymec_free(key);
  point_obs_t* data = val;
  polymec_free(data);
}

void model_define_point_observation(model_t* model, 
                                    const char* name, 
                                    real_t (*compute_point_observation)(void* context, 
                                                                        point_t* x,
                                                                        real_t t),
                                    point_t* point)

{
  point_obs_t* data = polymec_malloc(sizeof(point_obs_t));
  data->func = compute_point_observation;
  data->point = *point;
  string_ptr_unordered_map_insert_with_kv_dtor(model->point_obs, 
                                               string_dup(name), 
                                               data, point_obs_dtor);
}

typedef struct 
{
  real_t (*func)(void*, real_t);
} global_obs_t;

static void global_obs_dtor(char* key, void* val)
{
  polymec_free(key);
  polymec_free(val);
}

void model_define_global_observation(model_t* model, 
                                     const char* name, 
                                     real_t (*compute_global_observation)(void* context,
                                                                          real_t t))
{
  // ISO C doesn't allow us to stuff a function pointer into a pointer slot, so we
  // create a container for the function pointer and stuff it in there.
  // "Every problem can be solved by a layer of indirection," etc...
  global_obs_t* data = polymec_malloc(sizeof(global_obs_t));
  data->func = compute_global_observation;
  string_ptr_unordered_map_insert_with_kv_dtor(model->global_obs, 
                                               string_dup(name), 
                                               data, global_obs_dtor);
}

// Support for built-in observations.
static bool is_built_in_observation(const char* observation)
{
  return (!strcmp("sim_speed", observation) || 
          !strcmp("dt", observation) || 
          !strcmp("wall_time", observation));
}

static real_t built_in_observation(model_t* model, char* observation)
{
  if (!strcmp("sim_speed", observation))
    return model->sim_speed;
  else if (!strcmp("dt", observation))
    return model->dt;
  else if (!strcmp("wall_time", observation))
    return model->wall_time;
  else
  {
    polymec_error("Invalid built-in observation: %s", observation);
    return 0.0;
  }
}

static void key_dtor(char* key)
{
  string_free(key);
}

void model_observe(model_t* model, const char* observation)
{
  // Make sure this is an observation in the model.
  if (!is_built_in_observation(observation) &&
      !string_ptr_unordered_map_contains(model->point_obs, (char*)observation) && 
      !string_ptr_unordered_map_contains(model->global_obs, (char*)observation))
  {
    polymec_error("Observation '%s' is not defined for this model.", observation);
  }
  string_array_append_with_dtor(model->observations, string_dup(observation), key_dtor);
}

void model_set_observation_times(model_t* model, real_t* times, int num_times)
{
  ASSERT(num_times > 0);

  if (model->obs_times != NULL)
    polymec_free(model->obs_times);
  model->num_obs_times = num_times;
  model->obs_times = polymec_malloc(sizeof(real_t) * num_times);
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

  // Now record any observations we need to, given that the time step makes 
  // allowances for observations.
  if (model->num_obs_times > 0)
  {
    int obs_time_index = real_lower_bound(model->obs_times, model->num_obs_times, model->time);
    if ((obs_time_index < model->num_obs_times) && 
        (fabs(model->time - model->obs_times[obs_time_index]) < 1e-12)) // FIXME: Good enough?
      model_record_observations(model);
  }
}

// Initialize the model at the given time.
void model_init(model_t* model, real_t t)
{
  log_detail("%s: Initializing at time %g.", model->name, t);
  model->vtable.init(model->context, t);
  model->step = 0;
  model->dt = 0.0;
  model->time = t;
  model->wall_time0 = MPI_Wtime();
  model->wall_time = MPI_Wtime();
}

// Returns the largest permissible time step that can be taken by the model
// starting at time t.
real_t model_max_dt(model_t* model, char* reason)
{
  real_t dt = FLT_MAX;
  strcpy(reason, "No time step constraints.");

  // If we specified a maximum timestep, apply it here.
  if (model->max_dt < FLT_MAX)
  {
    dt = model->max_dt;
    strcpy(reason, "Max dt set in options.");
  }

  // If we have an observation time coming up, perhaps the next one will 
  // constrain the timestep.
  int obs_time_index = real_lower_bound(model->obs_times, model->num_obs_times, model->time);
  if (obs_time_index < model->num_obs_times)
  {
    real_t obs_time = model->obs_times[obs_time_index];
    real_t obs_dt = obs_time - model->time;
    ASSERT(obs_dt >= 0.0);
    if (obs_dt < 1e-12) 
    {
      // We're already at one observation time; set our sights on the next.
      if (obs_time_index < (model->num_obs_times - 1))
      {
        obs_time = model->obs_times[obs_time_index+1];
        dt = obs_time - model->time;
        sprintf(reason, "Requested observation time: %g", obs_time);
      }
    }
    else if (obs_dt < model->max_dt)
    {
      dt = obs_dt;
      sprintf(reason, "Requested observation time: %g", obs_time);
    }
    else if (2.0 * obs_dt < model->max_dt)
    {
      dt = obs_dt;
      sprintf(reason, "Requested observation time: %g", obs_time);
    }
  }

  // Now let the model have at it.
  if (model->vtable.max_dt != NULL)
    dt = model->vtable.max_dt(model->context, model->time, reason);
  return dt;
}

real_t model_advance(model_t* model, real_t max_dt)
{
  real_t pre_wall_time = MPI_Wtime();
  model->dt = model->vtable.advance(model->context, max_dt, model->time);

  // Check the time step.
  ASSERT(model->dt > 0.0);
  ASSERT(model->dt <= max_dt);

  model->time += model->dt;
  log_info("%s: Step %d (t = %g, dt = %g)", model->name, model->step, model->time, model->dt);
  model->step += 1;

  // Perform any busywork.
  model_do_periodic_work(model);

  // Record timings.
  real_t post_wall_time = MPI_Wtime();
  if (post_wall_time > pre_wall_time)
    model->sim_speed = model->dt / (post_wall_time - pre_wall_time); // Simulation "speed"
  else
    model->sim_speed = FLT_MAX;
  model->wall_time = post_wall_time;

  return model->dt;
}

real_t model_time(model_t* model)
{
  return model->time;
}

void model_finalize(model_t* model)
{
  log_detail("%s: Finalizing model at t = %g", model->name, model->time);
  if (model->vtable.finalize != NULL)
    model->vtable.finalize(model->context, model->step, model->time);
}

void model_load(model_t* model, int step)
{
  int nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  ASSERT(step >= 0);
  if (model->sim_name == NULL)
    polymec_error("No simulation name was set with model_set_sim_name.");
  char prefix[FILENAME_MAX], dir[FILENAME_MAX];
  snprintf(prefix, FILENAME_MAX, "%s", model->sim_name);
  snprintf(dir, FILENAME_MAX, "%s-%d", model->sim_name, nprocs);
  log_detail("%s: Loading save file from directory %s...", model->name, dir);
  model->vtable.load(model->context, prefix, dir, &model->time, step);
  model->step = step;

  // Reset the wall time(s).
  model->wall_time0 = MPI_Wtime();
  model->wall_time = MPI_Wtime();
}

void model_save(model_t* model)
{
  int nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  if (model->sim_name == NULL)
    polymec_error("No simulation name was set with model_set_sim_name.");
  char prefix[FILENAME_MAX], dir[FILENAME_MAX];
  snprintf(prefix, FILENAME_MAX, "%s", model->sim_name);
  snprintf(dir, FILENAME_MAX, "%s-%d", model->sim_name, nprocs);
  log_detail("%s: Writing save file to directory %s...", model->name, dir);
  model->vtable.save(model->context, prefix, dir, model->time, model->step);
}

void model_plot(model_t* model)
{
  int nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  if (model->sim_name == NULL)
    polymec_error("No simulation name was set with model_set_sim_name.");
  char prefix[FILENAME_MAX], dir[FILENAME_MAX];
  snprintf(prefix, FILENAME_MAX, "%s", model->sim_name);
  snprintf(dir, FILENAME_MAX, "%s-%d", model->sim_name, nprocs);
  log_detail("%s: Writing plot to directory %s...", model->name, dir);
  model->vtable.plot(model->context, prefix, dir, model->time, model->step);
}

void model_record_observations(model_t* model)
{
  if (model->sim_name == NULL)
    polymec_error("No simulation name was set with model_set_sim_name.");

  if ((model->observations->size == 0) || (model->num_obs_times == 0))
    return;

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Open up the observation file.
  char obs_fn[strlen(model->sim_name) + 5];
  snprintf(obs_fn, strlen(model->sim_name) + 5, "%s.obs", model->sim_name);
  FILE* obs = (rank == 0) ? fopen(obs_fn, "a") : NULL;

  // If we haven't recorded any observations yet, write a header.
  if ((model->time <= model->obs_times[0]) && (rank == 0))
  {
    log_detail("%s: Opening observation file %s.obs and writing header...", model->name, model->sim_name);

    // Provenance-related header.
    fprintf(obs, "# %s\n", polymec_invocation());
    time_t invoc_time = polymec_invocation_time();
    fprintf(obs, "# Invoked on: %s\n", ctime(&invoc_time));

    // Observation quantities.
    fprintf(obs, "# Observations for %s\n", model->sim_name);
    fprintf(obs, "# time ");
    for (int i = 0; i < model->observations->size; ++i)
      fprintf(obs, "%s ", model->observations->data[i]);
    fprintf(obs, "\n");
  }

  log_detail("%s: Recording observations at t = %g...", model->name, model->time);

  // Write the current simulation time.
  if (rank == 0)
    fprintf(obs, "%g ", model->time);

  // Go through all the desired observations.
  for (int i = 0; i < model->observations->size; ++i)
  {
    char* obs_name = model->observations->data[i];

    // Is this a built-in observation?
    real_t value;
    if (is_built_in_observation(obs_name))
    {
      value = built_in_observation(model, obs_name);
    }
    else
    {
      // Is it a point observation?
      point_obs_t** point_obs_data = (point_obs_t**)string_ptr_unordered_map_get(model->point_obs, obs_name);
      if (point_obs_data != NULL)
        value = (*point_obs_data)->func(model->context, &((*point_obs_data)->point), model->time);
      else
      {
        // It must be a global observation.
        global_obs_t** global_obs_data = (global_obs_t**)string_ptr_unordered_map_get(model->global_obs, obs_name);
        ASSERT(global_obs_data != NULL);
        value = (*global_obs_data)->func(model->context, model->time);
      }
    }

    // Write the observation to the file.
    if (rank == 0)
      fprintf(obs, "%g ", value);
  }
  if (rank == 0)
    fprintf(obs, "\n");

  // Close the file.
  if (rank == 0)
    fclose(obs);
}

void model_compute_error_norms(model_t* model, st_func_t* solution, real_t* error_norms)
{ 
  if (model->vtable.compute_error_norms != NULL)
    model->vtable.compute_error_norms(model->context, solution, model->time, error_norms);
  else
    polymec_error("%s: This model is not equipped to compute error norms.", model->name);
}

// This helper overrides interpreted parameters with options from the 
// command line.
static void override_interpreted_values(model_t* model, 
                                        real_t* t1, 
                                        real_t* t2, 
                                        int* max_steps)
{
  options_t* options = options_argv();
  // Run parameters -- not intrinsically part of a model.
  char* opt = options_value(options, "t1");
  if (opt != NULL)
    *t1 = atof(opt);
  opt = options_value(options, "t2");
  if (opt != NULL)
    *t2 = atof(opt);
  opt = options_value(options, "max_steps");
  if (opt != NULL)
    *max_steps = atoi(opt);

  opt = options_value(options, "load_step");
  if (opt != NULL)
  {
    model->load_step = atoi(opt);
    if (model->load_step < 0)
      polymec_error("Invalid load step: %d\n", model->load_step);
  }

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
  
  // Plot interval.
  char* plot_every = options_value(options, "plot_every");
  if (plot_every != NULL)
  {
    model->plot_every = atoi(plot_every);
    if (model->plot_every < 1)
      polymec_error("Invalid (non-positive) plot interval: %d\n", model->plot_every);
  }

  // Save interval.
  char* save_every = options_value(options, "save_every");
  if (save_every != NULL)
  {
    model->save_every = atoi(save_every);
    if (model->save_every < 1)
      polymec_error("Invalid (non-positive) save interval: %d\n", model->save_every);
  }

  // Handle observation times. Times can be specified with a regular interval
  // or a comma-delimited list.
  char* observe_every = options_value(options, "observe_every");
  if (observe_every != NULL)
  {
    model->observe_every = (real_t)atof(observe_every);
    if (model->observe_every <= 0.0)
      polymec_error("Invalid (non-positive) observation interval: %g\n", model->observe_every);
  }
  char* obs_times_str = options_value(options, "observation_times");
  if (obs_times_str != NULL)
  {
    if (model->observe_every > 0.0)
      polymec_error("Only one of observe_every and observation_times may be specified.");
    real_t* obs_times;
    int num_obs_times;
    obs_times = parse_observation_times(obs_times_str, &num_obs_times);
    if (obs_times != NULL)
    {
      model_set_observation_times(model, obs_times, num_obs_times);
      polymec_free(obs_times);
    }
    else
      polymec_error("Could not parse observation times string: %s\n", obs_times_str);
  }

  // If observation names are given, handle them here.
  string_array_clear(model->observations);
  char* obs_names_str = options_value(options, "observations");
  if (obs_names_str != NULL)
  {
    int num_obs;
    char** obs_names = string_split(obs_names_str, ",", &num_obs);
    for (int i = 0; i < num_obs; ++i)
    {
      model_observe(model, (const char*)obs_names[i]);
      polymec_free(obs_names[i]);
    }
    polymec_free(obs_names); 
  }

  char* max_dt = options_value(options, "max_dt");
  if (max_dt != NULL)
  {
    model->max_dt = atof(max_dt);
    if (model->max_dt <= 0.0)
      polymec_error("Invalid value for max_dt: %g", model->max_dt);
  }

  char* sim_name = options_value(options, "sim_name");
  if (sim_name != NULL)
    model_set_sim_name(model, sim_name);
}

void model_run(model_t* model, real_t t1, real_t t2, int max_steps)
{
  // Override options given at the command line.
  override_interpreted_values(model, &t1, &t2, &max_steps);

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
    // If we've not set observation times at this point, and if we got a 
    // value for observe_every, set the observation times.
    if (model->observe_every > 0.0)
    {
      int num_obs_times = (t2 - t1) / model->observe_every + 1;
      real_t* obs_times = polymec_malloc(sizeof(real_t) * num_obs_times);
      for (int i = 0; i < num_obs_times; ++i)
        obs_times[i] = i * model->observe_every;
      model_set_observation_times(model, obs_times, num_obs_times);
      polymec_free(obs_times);

      // Kick off the first set of observations!
      model_do_periodic_work(model);
    }
    else if (model->plot_every > 0)
      model_do_periodic_work(model);

    // Now run the calculation.
    while ((model->time < t2) && (model->step < max_steps))
    {
      char reason[POLYMEC_MODEL_MAXDT_REASON_SIZE];
      real_t max_dt = model_max_dt(model, reason);
      if (max_dt > t2 - model->time)
      {
        max_dt = t2 - model->time;
        snprintf(reason, POLYMEC_MODEL_MAXDT_REASON_SIZE, "End of simulation");
      }
      log_detail("%s: Max time step max_dt = %g\n (Reason: %s)", model->name, max_dt, reason);
      model_advance(model, max_dt);
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

static void driver_usage(const char* model_name)
{
  print_to_rank0("%s: usage:\n", model_name);
  print_to_rank0("%s [command] [args]\n\n", model_name);
  print_to_rank0("Here, [command] [args] is one of the following:\n\n");
  print_to_rank0("  run [file]                -- Runs a simulation with the given input file.\n");
  print_to_rank0("  benchmark [name]          -- Runs or queries a benchmark problem.\n");
  print_to_rank0("  help                      -- Prints information about the given model.\n\n");
  print_to_rank0("Benchmark commands:\n");
  print_to_rank0("  benchmark list            -- Lists all available benchmark problems.\n");
  print_to_rank0("  benchmark describe [name] -- Describes the given benchmark problem.\n");
  print_to_rank0("  benchmark all             -- Runs all available benchmark problems.\n");
  print_to_rank0("  benchmark [name]          -- Runs the given benchmark problem.\n\n");
  print_to_rank0("Help commands:\n");
  print_to_rank0("  help                      -- Prints model-specific help information.\n");
  print_to_rank0("  help list                 -- Prints a list of available functions.\n");
  print_to_rank0("  help [function/symbol]    -- Prints documentation for a function/symbol.\n");
  exit(-1);
}

// Prints model-specific help.
static void model_help(const char* exe_name, model_t* model, const char* arg)
{
  // If no argument was given, just print the model's basic documentation.
  if (arg == NULL)
  {
    if (model->doc != NULL)
    {
      int pos = 0;
      char* line;
      while (docstring_next(model->doc, &pos, &line))
        print_to_rank0("%s\n", line);
    }
    else
    {
      print_to_rank0("No documentation is available for the %s model.\n", model_name(model));
      print_to_rank0("Use '%s help list' to list available functions, and \n", model_name(model));
      print_to_rank0("'%s help <function>' for documentation on a given function.\n", model_name(model));
    }
  }
  else if (model->vtable.read_input != NULL)
  {
    static int rank = -1;
    if (rank == -1)
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0)
    {
      // Attempt to dig up the documentation for the given registered function.
      interpreter_t* interp = model_interpreter(model);
      interpreter_help(interp, arg, stderr);
    }
  }
  else
  {
    print_to_rank0("No specific help is available for models with custom input.\n");
  }
}

void model_set_sim_name(model_t* model, const char* sim_name)
{
  ASSERT(sim_name != NULL);
  if (model->sim_name != NULL)
    polymec_free(model->sim_name);
  model->sim_name = string_dup(sim_name);
}

int model_main(const char* model_name, model_ctor constructor, int argc, char* argv[])
{
  // Start everything up.
  polymec_init(argc, argv);

  // Get the parsed command line options.
  options_t* opts = options_argv();

  // Extract the executable name.
  char* full_exe_path = options_argument(opts, 0);
  char exe_dir[FILENAME_MAX], exe_name[FILENAME_MAX];
  parse_path(full_exe_path, exe_dir, exe_name);

  // Extract the command and arguments.
  char* command = options_argument(opts, 1);
  char* input = options_argument(opts, 2);

  // Validate our inputs.
  if (command != NULL)
  {
    int c = 0;
    static const char* valid_commands[] = {"run", "benchmark", "help", NULL};
    while (valid_commands[c] != NULL)
    {
      if (!strcmp(command, valid_commands[c]))
        break;
      ++c;
    }
    if (valid_commands[c] == NULL)
    {
      print_to_rank0("%s: invalid command: '%s'\n", exe_name, command);
      return -1;
    }
  }
  else
    driver_usage(exe_name);

  // Attempt to construct the model.
  model_t* model = (*constructor)();
  ASSERT(model != NULL);

  // Have we been asked for help?
  if (!strcmp(command, "help"))
  {
    model_help(exe_name, model, input);
    model_free(model);
    return 0;
  }

  // Have we been asked to do something related to a benchmark?
  if (!strcmp(command, "benchmark"))
  {
    if (input == NULL)
    {
      print_to_rank0("%s: No benchmark problem given! Usage:\n", exe_name);
      print_to_rank0("%s benchmark [problem] OR \n", exe_name);
      print_to_rank0("%s benchmark all OR \n", exe_name);
      print_to_rank0("%s benchmark list OR \n", exe_name);
      print_to_rank0("%s benchmark describe [problem] OR \n", exe_name);
      return -1;
    }

    if (!strcmp(input, "all")) // Run all benchmarks?
      model_run_all_benchmarks(model);
    else if (!strcmp(input, "list")) // List benchmarks?
    {
      print_to_rank0("Benchmarks for %s model:\n", model_name);
      int pos = 0;
      char *benchmark;
      model_benchmark_t* metadata;
      while (model_benchmark_map_next(model->benchmarks, &pos, &benchmark, &metadata))
      {
        if (metadata->description != NULL)
          print_to_rank0("  %s (%s)\n", benchmark, docstring_first_line(metadata->description));
        else
          print_to_rank0("  %s\n", benchmark);
      }
      print_to_rank0("\n");
    }
    else if (!strcmp(input, "describe")) // Describe a benchmark?
    {
      // We need to reach into argv for the benchmark name.
      if (argc < 4)
        print_to_rank0("%s: No benchmark specified for description.", exe_name);
      static int rank = -1;
      if (rank == -1)
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      if (rank == 0)
      {
        const char* benchmark = argv[3];
        model_describe_benchmark(model, benchmark, stderr);
      }
    }
    else
      model_run_benchmark(model, input);

    model_free(model);
    return 0;
  }

  // We are asked to run a simulation.
  ASSERT(!strcmp(command, "run"));
  if (input == NULL)
  {
    print_to_rank0("%s: No input file given! Usage:\n", exe_name);
    print_to_rank0("%s run [input file]\n", exe_name);
    return -1;
  }

  // Check to see whether the given file exists.
  FILE* fp = fopen(input, "r");
  if (fp == NULL)
  {
    print_to_rank0("%s: Input file not found: %s\n", exe_name, input);
    return -1;
  }
  fclose(fp);

  // By default, the simulation is named after the input file (minus its path 
  // and anything after the first alphanumeric character).
  char dir_name[FILENAME_MAX], file_name[FILENAME_MAX];
  parse_path(input, dir_name, file_name);
  int len = strlen(file_name), end = 0;
  while (end < len)
  {
    if (!isalnum(file_name[end]))
      break;
    ++end;
  }
  char prefix[FILENAME_MAX];
  strncpy(prefix, file_name, end);
  prefix[end] = '\0';
  model_set_sim_name(model, prefix);

  // Read the contents of the input file into the model's interpreter.
  model_read_input_file(model, input);

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

  // Run the model.
  model_run(model, t1, t2, max_steps);

  // Clean up.
  model_free(model);

  return 0;
}

static void multi_model_usage(const char* exe_name, 
                              model_dispatch_t models[])
{
  print_to_rank0("%s: usage:\n", exe_name);
  print_to_rank0("%s [command] [model] [args]\n\n", exe_name);
  print_to_rank0("Here, [command] [model] [args] is one of the following:\n\n");
  print_to_rank0("  run [model] [file]                -- Runs a simulation with the given\n");
  print_to_rank0("                                       input file.\n");
  print_to_rank0("  benchmark [model] [name]          -- Runs or queries a benchmark problem.\n");
  print_to_rank0("  help [model]                      -- Prints information about the\n");
  print_to_rank0("                                       given model.\n\n");
  print_to_rank0("Benchmark commands:\n");
  print_to_rank0("  benchmark [model] list            -- Lists all available benchmark problems.\n");
  print_to_rank0("  benchmark [model] describe [name] -- Describes the given benchmark problem.\n");
  print_to_rank0("  benchmark [model] all             -- Runs all available benchmark problems.\n");
  print_to_rank0("  benchmark [model] [name]          -- Runs the given benchmark problem.\n\n");
  print_to_rank0("Help commands:\n");
  print_to_rank0("  help [model]                      -- Prints model-specific help information.\n");
  print_to_rank0("  help [model] list                 -- Prints a list of available functions.\n");
  print_to_rank0("  help [model] [function/symbol]    -- Prints documentation for a\n");
  print_to_rank0("                                       function/symbol.\n\n");
  print_to_rank0("Above, [model] is one of:\n");
  int i = 0;
  while (strcmp(models[i].model_name, END_OF_MODELS.model_name) != 0)
  {
    print_to_rank0("  %s\n", models[i].model_name);
    ++i;
  }
  exit(-1);
}

// Prints model-specific help for multi-model programs.
static void multi_model_help(const char* exe_name, model_t* model, const char* arg)
{
  // If no argument was given, just print the model's basic documentation.
  if (arg == NULL)
  {
    if (model->doc != NULL)
    {
      int pos = 0;
      char* line;
      while (docstring_next(model->doc, &pos, &line))
        print_to_rank0("%s\n", line);
    }
    else
    {
      print_to_rank0("No documentation is available for the %s model.\n", model_name(model));
      print_to_rank0("Use '%s help %s list' to list available functions, and \n", exe_name, model_name(model));
      print_to_rank0("'%s help %s <function>' for documentation on a given function.\n", exe_name, model_name(model));
    }
  }
  else if (model->vtable.read_input != NULL)
  {
    // Attempt to dig up the documentation for the given registered function.
    static int rank = -1;
    if (rank == -1)
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0)
    {
      interpreter_t* interp = model_interpreter(model);
      interpreter_help(interp, arg, stderr);
    }
  }
  else
  {
    print_to_rank0("No specific help is available for models with custom input.\n");
  }
}


int multi_model_main(model_dispatch_t model_table[], 
                     int argc, 
                     char* argv[])
{
  // Start everything up.
  polymec_init(argc, argv);

  // Get the parsed command line options.
  options_t* opts = options_argv();

  // Extract the executable name.
  char* full_exe_path = options_argument(opts, 0);
  char exe_dir[FILENAME_MAX], exe_name[FILENAME_MAX];
  parse_path(full_exe_path, exe_dir, exe_name);

  // Extract the command and arguments.
  char* command = options_argument(opts, 1);
  char* model_name = options_argument(opts, 2);
  char* input = options_argument(opts, 3);

  // Validate our inputs.
  if (command != NULL)
  {
    int c = 0;
    static const char* valid_commands[] = {"run", "benchmark", "help", NULL};
    while (valid_commands[c] != NULL)
    {
      if (!strcmp(command, valid_commands[c]))
        break;
      ++c;
    }
    if (valid_commands[c] == NULL)
    {
      print_to_rank0("%s: invalid command: '%s'\n", exe_name, command);
      return -1;
    }

    // Show generic help for help with no model name.
    if ((strcmp(command, "help") == 0) && (model_name == NULL))
      multi_model_usage(exe_name, model_table);
  }
  else
    multi_model_usage(exe_name, model_table);
  if (model_name == NULL)
    multi_model_usage(exe_name, model_table);

  // Attempt to construct the model.
  model_t* model = NULL;
  {
    int i = 0;
    while (strcmp(model_table[i].model_name, END_OF_MODELS.model_name) != 0)
    {
      if (strcmp(model_table[i].model_name, model_name) == 0)
      {
        model = model_table[i].model_constructor();
        break;
      }
      ++i;
    }
  }
  if (model == NULL)
  {
    print_to_rank0("%s: Invalid model: %s\n", exe_name, model_name);
    print_to_rank0("%s: Valid models are:\n", exe_name);
    int i = 0;
    while (strcmp(model_table[i].model_name, END_OF_MODELS.model_name) != 0)
    {
      print_to_rank0("  %s\n", model_table[i].model_name);
      ++i;
    }
    return -1;
  }

  // Have we been asked for help?
  if (!strcmp(command, "help"))
  {
    multi_model_help(exe_name, model, input);
    model_free(model);
    return 0;
  }

  // Have we been asked to do something related to a benchmark?
  if (!strcmp(command, "benchmark"))
  {
    if (input == NULL)
    {
      print_to_rank0("%s: No benchmark problem given! Usage:\n", exe_name);
      print_to_rank0("%s benchmark [model] [problem] OR \n", exe_name);
      print_to_rank0("%s benchmark [model] all OR \n", exe_name);
      print_to_rank0("%s benchmark [model] list OR \n", exe_name);
      print_to_rank0("%s benchmark [model] describe [problem] OR \n", exe_name);
      return -1;
    }

    if (!strcmp(input, "all")) // Run all benchmarks?
      model_run_all_benchmarks(model);
    else if (!strcmp(input, "list")) // List benchmarks?
    {
      print_to_rank0("Benchmarks for %s model:\n", model_name);
      int pos = 0;
      char *benchmark;
      model_benchmark_t* metadata;
      while (model_benchmark_map_next(model->benchmarks, &pos, &benchmark, &metadata))
      {
        if (metadata->description != NULL)
          print_to_rank0("  %s (%s)\n", benchmark, docstring_first_line(metadata->description));
        else
          print_to_rank0("  %s\n", benchmark);
      }
      print_to_rank0("\n");
    }
    else if (!strcmp(input, "describe")) // Describe a benchmark?
    {
      // We need to reach into argv for the benchmark name.
      if (argc < 4)
        print_to_rank0("%s: No benchmark specified for description.", exe_name);
      static int rank = -1;
      if (rank == -1)
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      if (rank == 0)
      {
        const char* benchmark = argv[3];
        model_describe_benchmark(model, benchmark, stderr);
      }
    }
    else
      model_run_benchmark(model, input);

    model_free(model);
    return 0;
  }

  // We are asked to run a simulation.
  ASSERT(!strcmp(command, "run"));
  if (input == NULL)
  {
    print_to_rank0("%s: No input file given! Usage:\n", exe_name);
    print_to_rank0("%s run %s [input file]\n", exe_name, model_name);
    return -1;
  }

  // Check to see whether the given file exists.
  FILE* fp = fopen(input, "r");
  if (fp == NULL)
  {
    print_to_rank0("%s: Input file not found: %s\n", exe_name, input);
    return -1;
  }
  fclose(fp);

  // By default, the simulation is named after the input file (minus its path 
  // and anything after the first alphanumeric character).
  char dir_name[FILENAME_MAX], file_name[FILENAME_MAX];
  parse_path(input, dir_name, file_name);
  int len = strlen(file_name), end = 0;
  while (end < len)
  {
    if (!isalnum(file_name[end]))
      break;
    ++end;
  }
  char prefix[FILENAME_MAX];
  strncpy(prefix, file_name, end);
  prefix[end] = '\0';
  model_set_sim_name(model, prefix);

  // Read the contents of the input file into the model's interpreter.
  model_read_input_file(model, input);

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

  // Run the model.
  model_run(model, t1, t2, max_steps);

  // Clean up.
  model_free(model);

  return 0;
}

static void minimal_driver_usage(const char* model_name)
{
  print_to_rank0("%s: usage:\n", model_name);
  print_to_rank0("%s input [option1=val [option2=val2 [...]]]\n\n", model_name);
  exit(-1);
}

int model_minimal_main(const char* model_name, model_ctor constructor, int argc, char* argv[])
{
  // Start everything up.
  polymec_init(argc, argv);

  // Get the parsed command line options.
  options_t* opts = options_argv();

  // Here, the first argument serves as the input.
  char* input = options_argument(opts, 1);
  if ((input == NULL) || !strcmp(input, "help"))
  {
    minimal_driver_usage(model_name);
  }

  // Attempt to construct the model.
  model_t* model = (*constructor)();
  ASSERT(model != NULL);

  // Check to see whether the given file exists.
  FILE* fp = fopen(input, "r");
  if (fp == NULL)
  {
    print_to_rank0("%s: Input file not found: %s\n", model_name, input);
    return -1;
  }
  fclose(fp);

  // By default, the simulation is named after the input file (minus its suffix).
  model_set_sim_name(model, input);

  // Read the contents of the input file into the model's interpreter.
  model_read_input_file(model, input);

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

  // Run the model.
  model_run(model, t1, t2, max_steps);

  // Clean up.
  model_free(model);

  return 0;
}

void model_report_conv_rate(real_t conv_rate, real_t sigma)
{
  options_t* options = options_argv();
  char crstr[1024], sigmastr[1024];
  snprintf(crstr, 1024, "%g", conv_rate);
  snprintf(sigmastr, 1024, "%g", sigma);
  options_set(options, "conv_rate", crstr);
  options_set(options, "conv_rate_sigma", sigmastr);
}

void model_report_error_norm(real_t error_norm)
{
  options_t* options = options_argv();
  char nstr[1024];
  snprintf(nstr, 1024, "%g", error_norm);
  options_set(options, "error_norm", nstr);
}

