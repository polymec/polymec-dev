// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <stdarg.h>
#include <signal.h>
#include "core/polymec.h"
#include "core/unordered_set.h"
#include "core/unordered_map.h"
#include "core/array.h"
#include "core/array_utils.h"
#include "core/text_buffer.h"
#include "core/timer.h"
#include "core/hash_functions.h"
#include "core/lua_driver.h"
#include "model/model.h"

#ifdef _OPENMP
#include <omp.h>
#endif

static int _mpi_rank = -1;

// This type maps a probe to an array of acquisition times.
static inline int probe_hash(probe_t* p)
{
  return string_hash(probe_name(p));
}

static inline bool probe_equals(probe_t* p, probe_t* q)
{
  return ptr_equals(p, q);
}

DEFINE_UNORDERED_MAP(probe_map, probe_t*, real_array_t*, probe_hash, probe_equals)
DEFINE_UNORDERED_MAP(probe_data_map, char*, probe_data_array_t*, string_hash, string_equals)

// This catches the SIGINT (Ctrl-C) signal and sets a flag for the model
// to respond appropriately.
typedef void (*sighandler_t)(int sig);
static int _received_signal = 0;
static bool _signal_handlers_set = false;
static void handle_sigint_or_sigterm(int signal)
{
  _received_signal = signal;
}

struct model_t
{
  // Model metadata.
  void* context;
  char* name;
  model_vtable vtable;
  model_parallelism_t parallelism;

  int save_every;       // Save frequency (steps).
  real_t plot_every;    // Plot frequency (time units).
  bool plot_this_step;  // Flag to plot after the current step completes.

  int load_step; // -1 if starting a new sim, >= 0 if loading from a file.

  // Probes and their acquisition times.
  probe_map_t* probes;

  // Acquired probe data.
  probe_data_map_t* probe_data;

  // Diagnostics mode.
  model_diag_mode_t diag_mode;

  // Intercept SIGINT and SIGTERM?
  bool handle_signals;

  // Data related to a given simulation.
  char* sim_prefix; // Simulation naming prefix.
  char* sim_dir;    // Simulation directory.
  real_t time;       // Current simulation time.
  real_t wall_time0; // Wall time at simulation start.
  real_t wall_time;  // Current wall time.
  real_t sim_speed;  // Simulation "speed" (sim time / wall time per step).
  int step;          // Current simulation step number.
  real_t dt;         // Current simulation time step.
  real_t initial_dt; // Initial time step.
  real_t max_dt;     // Maximum time step.
  real_t min_dt;     // Minimum time step.
};

// Here's a static set of all named model instances.
static string_unordered_set_t* model_singletons = NULL;

// This helper ensures that the given model is the only instance of its
// type in memory if it has parallelism type MPI_SERIAL_SINGLETON or
// MPI_PARALLEL_SINGLETON.
static void enforce_singleton_instances(model_t* model)
{
  if ((model->parallelism == MODEL_SERIAL_SINGLETON) ||
      (model->parallelism == MODEL_MPI_SINGLETON))
  {
    if (model_singletons == NULL)
      model_singletons = string_unordered_set_new();
    if (string_unordered_set_contains(model_singletons, model->name))
      polymec_error("Model %s is a singleton: only 1 instance can exist in memory.", model->name);

    // NOTE: the model itself outlives its entry in the list of singletons,
    // NOTE: so we don't need to insert a copy of the name.
    string_unordered_set_insert(model_singletons, model->name);
  }
}

static void enforce_mpi_parallelism(model_t* model)
{
  if ((model->parallelism == MODEL_SERIAL) ||
      (model->parallelism == MODEL_SERIAL_SINGLETON))
  {
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    if (nprocs > 1)
      polymec_error("Model %s is a serial model and can only be run on a single process.", model->name);
  }
}

// This helper checks the given model for thread safety.
static void check_thread_safety(model_t* model)
{
#ifdef _OPENMP
  int num_threads = omp_get_num_threads();
#else
  int num_threads = 1;
#endif
  if ((num_threads > 1) && (model->parallelism != MODEL_MPI_THREAD_SAFE))
  {
    polymec_error("Model %s is not thread-safe but is being invoked with %d threads.",
                  model->name, num_threads);
  }
}

model_t* model_new(const char* name,
                   void* context,
                   model_vtable vtable,
                   model_parallelism_t parallelism)
{
  if (_mpi_rank == -1)
    MPI_Comm_rank(MPI_COMM_WORLD, &_mpi_rank);

  model_t* model = polymec_malloc(sizeof(model_t));
  model->vtable = vtable;
  model->context = context;
  model->name = string_dup(name);
  model->parallelism = parallelism;
  model->save_every = -1;
  model->plot_every = -REAL_MAX;
  model->load_step = -1;
  model->wall_time = 0.0;
  model->wall_time0 = 0.0;
  model->initial_dt = REAL_MAX;
  model->time = 0.0;
  model->dt = 0.0;
  model->step = 0;
  model->max_dt = REAL_MAX;
  model->min_dt = 0.0;
  model->diag_mode = MODEL_DIAG_NEAREST_STEP;

  // By default, we make model steps uninterruptible by intercepting
  // SIGINT and SIGTERM.
  model->handle_signals = true;

  // Set defaults for the simulation file prefix and directory.
  const char* script = lua_driver_script();
  if (script != NULL)
  {
    char dirname[FILENAME_MAX], filename[FILENAME_MAX];
    parse_path(script, dirname, filename);
    model->sim_prefix = string_dup(filename);
    model->sim_dir = string_dup(filename);
  }
  else
  {
    model->sim_prefix = string_dup(name);
    model->sim_dir = string_dup(".");
  }

  // Initialize probe and probe data maps.
  model->probes = probe_map_new();
  model->probe_data = probe_data_map_new();

  // Enforce our given parallelism model.
  enforce_singleton_instances(model);
  enforce_mpi_parallelism(model);

  return model;
}

void model_free(model_t* model)
{
  // If the model is a singleton, we remove its instance from
  // our set of running singletons.
  if ((model->parallelism == MODEL_SERIAL_SINGLETON) ||
      (model->parallelism == MODEL_MPI_SINGLETON))
  {
    string_unordered_set_delete(model_singletons, model->name);

    // If this is the last singleton, delete our set of instances.
    if (string_unordered_set_empty(model_singletons))
    {
      string_unordered_set_free(model_singletons);
      model_singletons = NULL;
    }
  }

  if ((model->context != NULL) && (model->vtable.dtor != NULL))
    model->vtable.dtor(model->context);
  polymec_free(model->name);

  polymec_free(model->sim_prefix);
  polymec_free(model->sim_dir);

  // Clear probe stuff.
  probe_map_free(model->probes);
  probe_data_map_free(model->probe_data);

  polymec_free(model);
}

char* model_name(model_t* model)
{
  return model->name;
}

model_parallelism_t model_parallelism(model_t* model)
{
  return model->parallelism;
}

void model_handle_signals(model_t* model, bool flag)
{
  model->handle_signals = flag;
}

extern void probe_set_model(probe_t* probe, model_t* model);

void model_add_probe(model_t* model,
                     probe_t* probe,
                     real_t* acq_times,
                     size_t num_acq_times)
{
  log_debug("%s: Adding probe %s (measures %s)", model->name,
            probe_name(probe), probe_data_name(probe));
  real_array_t* times = real_array_new();
  real_array_resize(times, num_acq_times);
  memcpy(times->data, acq_times, sizeof(real_t) * num_acq_times);
  probe_map_insert_with_kv_dtors(model->probes, probe, times, probe_free, real_array_free);

  // Perform any work that needs doing.
  if (model->vtable.add_probe != NULL)
    model->vtable.add_probe(model->context, probe_context(probe));
  probe_set_model(probe, model);
}

static void model_do_periodic_work(model_t* model)
{
  // Do plots and saves.
  if (model->plot_every > 0.0)
  {
    if (model->plot_this_step)
    {
      model_plot(model);
      model->plot_this_step = false;
    }
  }

  // Save if the step # is right and if we're not on a freshly-loaded step.
  if ((model->save_every > 0) &&
      ((model->step % model->save_every) == 0) &&
      (model->load_step != model->step))
    model_save(model);

  // Now acquire any data we need to, given that the time step makes
  // allowances for acquisitions.
  model_acquire(model);
}

// Initialize the model at the given time.
void model_init(model_t* model, real_t t)
{
  START_FUNCTION_TIMER();
  check_thread_safety(model);

  log_detail("%s: Initializing at time %g.", model->name, t);
  model->vtable.init(model->context, t);
  model->step = 0;
  model->dt = model->initial_dt;
  model->time = t;
  model->wall_time0 = MPI_Wtime();
  model->wall_time = MPI_Wtime();
  model->plot_this_step = true;
  STOP_FUNCTION_TIMER();
}

real_t model_initial_dt(model_t* model)
{
  return model->initial_dt;
}

void model_set_initial_dt(model_t* model, real_t dt0)
{
  ASSERT(dt0 > 0.0);
  model->initial_dt = dt0;
}

// Returns the next scheduled acquisition time.
static real_t model_next_acq_time(model_t* model)
{
  real_t acq_time = REAL_MAX;

  int pos = 0;
  probe_t* probe;
  real_array_t* acq_times;
  while (probe_map_next(model->probes, &pos, &probe, &acq_times))
  {
    size_t acq_time_index = real_lower_bound(acq_times->data, acq_times->size, model->time);
    if (acq_time_index < acq_times->size)
    {
      real_t t = acq_times->data[acq_time_index];
      real_t dt = t - model->time;
      ASSERT(dt >= 0.0);
      if (dt < 1e-12)
      {
        // We're already at one observation time; set our sights on the next.
        if (acq_time_index < (acq_times->size - 1))
        {
          acq_time = MIN(acq_time, acq_times->data[acq_time_index+1]);
          dt = t - model->time;
        }
      }
      else
        acq_time = MIN(acq_time, t);
    }
  }

  return acq_time;
}

real_t model_max_dt(model_t* model, char** reason)
{
  static char _reason[1025];
  real_t dt = REAL_MAX;
  strncpy(_reason, "No time step constraints.", 1024);

  // First let the model have at it.
  if (model->vtable.max_dt != NULL)
  {
    char model_reason[1025];
    real_t model_dt = model->vtable.max_dt(model->context, model->time, model_reason);
    if (model_dt < dt)
    {
      dt = model_dt;
      strncpy(_reason, model_reason, 1024);
    }
  }

  // If we specified a maximum timestep, apply it here.
  if (model->max_dt < dt)
  {
    dt = model->max_dt;
    strncpy(_reason, "Specified maximum timestep.", 1024);
  }

  // The following constraints only apply if we are collecting diagnostics
  // at exact times.
  if (model->diag_mode == MODEL_DIAG_EXACT_TIME)
  {
    // If we have an observation time coming up, perhaps the next one will
    // constrain the timestep.
    real_t acq_time = model_next_acq_time(model);
    real_t acq_dt = acq_time - model->time;
    if (acq_dt < dt)
    {
      dt = acq_dt;
      snprintf(_reason, 1024, "Requested acquisition time: %g", acq_time);
    }

    // If we have a plot time coming up, perhaps the next one will
    // constrain the timestep.
    if (model->plot_every > 0.0)
    {
      int n = (int)(model->time / model->plot_every);
      real_t plot_time = (n+1) * model->plot_every;
      real_t plot_dt = plot_time - model->time;
      ASSERT(plot_dt >= 0.0);
      if (plot_dt < 1e-12)
      {
        // We're already at one plot time; set our sights on the next.
        plot_dt = model->plot_every;
      }

      // Is our plot timestep the limiting factor?
      if (plot_dt < dt)
      {
        dt = plot_dt;
        model->plot_this_step = true;
        snprintf(_reason, 1024, "Requested plot time: %g", plot_time);
      }
      else if (2.0 * plot_dt < dt)
      {
        dt = plot_dt;
        model->plot_this_step = true;
        snprintf(_reason, 1024, "Requested plot time: %g", plot_time);
      }
      else if (reals_nearly_equal(plot_dt, dt, 1e-12)) // FIXME: corny
        model->plot_this_step = true;
    }
  }
  else
  {
    // We still need to create plots, even if we're not stopping exactly
    // at the right time.
    if (model->plot_every > 0.0)
    {
      int n = (int)(model->time / model->plot_every);
      real_t plot_time = (n+1) * model->plot_every;
      real_t plot_dt = plot_time - model->time;
      if (plot_dt < dt)
        model->plot_this_step = true;
    }
  }

  if (reason != NULL)
    *reason = _reason;

  return dt;
}

void model_set_max_dt(model_t* model, real_t max_dt)
{
  ASSERT(max_dt > 0.0);
  model->max_dt = max_dt;
}

void model_set_min_dt(model_t* model, real_t min_dt)
{
  ASSERT(min_dt >= 0.0);
  ASSERT(min_dt <= model->max_dt);
  model->min_dt = min_dt;
}

real_t model_min_dt(model_t* model)
{
  return model->min_dt;
}

real_t model_advance(model_t* model, real_t max_dt)
{
  START_FUNCTION_TIMER();

  // If we're not inside model_advance(), set up our SIGINT handler.
  sighandler_t prev_sigint_handler = NULL;
  sighandler_t prev_sigterm_handler = NULL;
  if (!_signal_handlers_set && model->handle_signals)
  {
    prev_sigint_handler = signal(SIGINT, handle_sigint_or_sigterm);
    prev_sigterm_handler = signal(SIGTERM, handle_sigint_or_sigterm);
  }
  else
  {
    // If we've received a SIGINT, exit immediately.
    if (_received_signal)
    {
      STOP_FUNCTION_TIMER();
      return 0.0;
    }
  }

  real_t pre_wall_time = MPI_Wtime();
  model->dt = model->vtable.advance(model->context, max_dt, model->time);

  // Check the time step.
  static real_t dt_fuzz = -1.0;
  if (dt_fuzz < 0.0)
    dt_fuzz = pow(REAL_EPSILON, 2.0/3.0);

  ASSERT(model->dt > 0.0);
  ASSERT(model->dt <= max_dt + dt_fuzz);

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
    model->sim_speed = REAL_MAX;
  model->wall_time = post_wall_time;

  // Restore the previous SIGINT handler if we set it.
  if (!_signal_handlers_set && model->handle_signals)
  {
    signal(SIGINT, prev_sigint_handler);
    signal(SIGTERM, prev_sigterm_handler);

    // If we got a signal, propagate it.
    if (_received_signal != 0)
    {
      _received_signal = 0;
      raise(_received_signal);
    }
  }

  STOP_FUNCTION_TIMER();
  return model->dt;
}

real_t model_time(model_t* model)
{
  return model->time;
}

int model_step(model_t* model)
{
  return model->step;
}

void model_finalize(model_t* model)
{
  START_FUNCTION_TIMER();
  log_detail("%s: Finalizing model at t = %g", model->name, model->time);
  if (model->vtable.finalize != NULL)
    model->vtable.finalize(model->context, model->step, model->time);

  // Do any postprocessing required by probes.
  int pos = 0;
  probe_t* probe;
  real_array_t* acq_times;
  while (probe_map_next(model->probes, &pos, &probe, &acq_times))
  {
    char* data_name = probe_data_name(probe);
    probe_data_array_t** data_p = probe_data_map_get(model->probe_data, data_name);
    probe_data_array_t* data = (data_p != NULL) ? *data_p : NULL;
    probe_postprocess(probe, acq_times, data);
  }

  // Clear any plotting, saving, loading that has been requested.
  model->plot_every = -REAL_MAX;
  model->load_step = -1;
  model->save_every = -1;

  STOP_FUNCTION_TIMER();
}

bool model_load(model_t* model, int step)
{
  START_FUNCTION_TIMER();
  check_thread_safety(model);

  // Is loading supported by this model?
  if (model->vtable.load == NULL)
    polymec_error("Loading from save files is not supported by this model.");

  ASSERT(step >= 0);
  log_detail("%s: Loading save file from directory %s...", model->name, model->sim_dir);
  bool loaded = model->vtable.load(model->context, model->sim_prefix, model->sim_dir, &model->time, step);
  if (loaded)
  {
    model->step = step;
    model->plot_this_step = true;

    // Reset the wall time(s).
    model->wall_time0 = MPI_Wtime();
    model->wall_time = MPI_Wtime();
  }
  else
    log_detail("%s: Could not load save file.", model->name);
  STOP_FUNCTION_TIMER();
  return loaded;
}

void model_save(model_t* model)
{
  START_FUNCTION_TIMER();
  // Is saving supported by this model?
  if (model->vtable.save == NULL)
    polymec_error("Saving is not supported by this model.");

  log_detail("%s: Writing save file to directory %s...", model->name, model->sim_dir);
  model->vtable.save(model->context, model->sim_prefix, model->sim_dir, model->time, model->step);
  STOP_FUNCTION_TIMER();
}

void model_plot(model_t* model)
{
  START_FUNCTION_TIMER();
  if (model->vtable.plot == NULL)
    polymec_error("Plotting is not supported by this model.");

  log_detail("%s: Writing plot to directory %s...", model->name, model->sim_dir);
  char plot_prefix[FILENAME_MAX+1];
  snprintf(plot_prefix, FILENAME_MAX, "%s-plot", model->sim_prefix);
  model->vtable.plot(model->context, plot_prefix, model->sim_dir, model->time, model->step);
  STOP_FUNCTION_TIMER();
}

void model_acquire(model_t* model)
{
  START_FUNCTION_TIMER();
  int pos = 0;
  probe_t* probe;
  real_array_t* acq_times;
  while (probe_map_next(model->probes, &pos, &probe, &acq_times))
  {
    // Figure out whether to actually acquire the data now.
    bool acquire_now = false;
    size_t acq_time_index = real_lower_bound(acq_times->data, acq_times->size,
                                             model->time);
    if ((acq_time_index < acq_times->size) &&
         reals_nearly_equal(model->time, acq_times->data[acq_time_index], 1e-12)) // FIXME: Good enough?
      acquire_now = true;
    else if (model->diag_mode == MODEL_DIAG_NEAREST_STEP)
    {
      real_t last_dt = model->dt;
      size_t last_acq_time_index = real_lower_bound(acq_times->data, acq_times->size,
                                                    model->time - last_dt);
      if (last_acq_time_index < acq_time_index)
        acquire_now = true;
    }

    // Do what needs doing.
    if (acquire_now)
    {
      // Acquire data from this probe.
      probe_data_t* data = probe_acquire(probe, model->time);

      if (_mpi_rank == 0)
      {
        // Get an array in which to stash this data.
        char* data_name = probe_data_name(probe);
        probe_data_array_t** array_p = probe_data_map_get(model->probe_data, data_name);
        probe_data_array_t* array = NULL;
        if (array_p == NULL)
        {
          array = probe_data_array_new();
          probe_data_map_insert_with_kv_dtors(model->probe_data,
                                              string_dup(data_name), array,
                                              string_free, probe_data_array_free);
        }
        else
          array = *array_p;

        // Stash it!
        probe_data_array_append_with_dtor(array, data, probe_data_free);
      }
      else
      {
        // Jettison the data on other ranks.
        if (data != NULL)
          probe_data_free(data);
      }
    }
  }
  STOP_FUNCTION_TIMER();
}

void model_run(model_t* model, real_t t1, real_t t2, int max_steps)
{
  START_FUNCTION_TIMER();
  ASSERT(t2 >= t1);

  sighandler_t prev_sigint_handler = NULL;
  sighandler_t prev_sigterm_handler = NULL;
  if (model->handle_signals)
  {
    // Set up our SIGINT handler.
    _received_signal = 0;
    prev_sigint_handler = signal(SIGINT, handle_sigint_or_sigterm);
    prev_sigterm_handler = signal(SIGTERM, handle_sigint_or_sigterm);
    _signal_handlers_set = true;
  }

  bool model_initialized = false;
  if (model->load_step == -1)
  {
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
    model_initialized = true;
  }
  else
  {
    model_initialized = model_load(model, model->load_step);
    t1 = model->time;
  }

  if (!model_initialized)
  {
    signal(SIGINT, prev_sigint_handler);
    signal(SIGTERM, prev_sigterm_handler);
    _signal_handlers_set = false;
    STOP_FUNCTION_TIMER();
    polymec_error("%s: Could not load model from step %d", model->name, model->load_step);
  }

  if (reals_equal(t1, t2))
    model_do_periodic_work(model);
  else
  {
    if ((model->plot_every > 0.0) || (model->save_every > 0))
      model_do_periodic_work(model);
    else
      model_acquire(model); // make sure probes can acquire data at t1.

    // Now run the calculation.
    while ((model->time < t2) && (model->step < max_steps) && (_received_signal == 0))
    {
      char reason[1025];

      // Let the model tell us the maximum time step it can take.
      char* model_reason;
      real_t max_dt = model_max_dt(model, &model_reason);
      strncpy(reason, model_reason, 1024);

      // Have we fallen below the minimum allowable time step? If so,
      // we end the simulation.
      if (max_dt < model->min_dt)
      {
        log_urgent("%s: Max time step (%g) has fallen below min time step (%g).", model->name, max_dt, model->min_dt);
        log_urgent("%s: Reason: %s", model->name, reason);
        log_urgent("%s: The simulation will now terminate.", model->name);
        break;
      }

      // If it's the first step, apply the initial time step size as a
      // constraint.
      if (model->step == 0)
      {
        if (model->initial_dt < max_dt)
        {
          max_dt = model->initial_dt;
          snprintf(reason, 1024, "Initial time step size");
        }
      }

      // If we are approaching the end of the simulation, we need to apply
      // the end time as a constraint.
      real_t remaining_time = t2 - model->time;
      if (max_dt > remaining_time)
      {
        max_dt = remaining_time;
        snprintf(reason, 1024, "End of simulation");
      }

      // If we are not quite there but we are close enough to encounter the
      // minimum time step, break the final 2 steps into more equally-spaced
      // sizes.
      else if ((max_dt + model->min_dt) > 2.0 * remaining_time)
      {
        max_dt = 0.5 * remaining_time;
        snprintf(reason, 1024, "End of simulation (subject to min_dt).");
      }

      log_detail("%s: Max time step max_dt = %g\n (Reason: %s)", model->name, max_dt, reason);
      model_advance(model, max_dt);
    }
    if (_received_signal != 0)
      log_urgent("%s: Simulation interrupted at step %d.", model->name, model->step);
  }

  if (model->handle_signals)
  {
    // Restore the previous SIGINT handler.
    signal(SIGINT, prev_sigint_handler);
    signal(SIGTERM, prev_sigterm_handler);
    _signal_handlers_set = false;

    // If we got a signal, propagate it.
    if (_received_signal != 0)
      raise(_received_signal);
  }

  // Do any finalization.
  model_finalize(model);
  log_detail("%s: Run concluded at time %g.", model->name, model->time);

  STOP_FUNCTION_TIMER();
}

void* model_context(model_t* model)
{
  return model->context;
}

bool model_next_probe_data(model_t* model,
                           int* pos,
                           char** quantity,
                           probe_data_array_t** data)
{
  return probe_data_map_next(model->probe_data, pos, quantity, data);
}

probe_data_array_t* model_probe_data(model_t* model, const char* quantity)
{
  probe_data_array_t** array_p = probe_data_map_get(model->probe_data, (char*)quantity);
  if (array_p != NULL)
    return *array_p;
  else
    return NULL;
}

void model_set_prefix(model_t* model, const char* prefix)
{
  ASSERT(prefix != NULL);
  log_debug("%s: Setting prefix to \"%s\".", model->name, prefix);
  polymec_free(model->sim_prefix);
  model->sim_prefix = string_dup(prefix);
}

void model_set_directory(model_t* model, const char* directory)
{
  ASSERT(directory != NULL);
  log_debug("%s: Setting directory to \"%s\".", model->name, directory);
  polymec_free(model->sim_dir);
  model->sim_dir = string_dup(directory);
}

void model_plot_every(model_t* model, real_t T)
{
  ASSERT(T > 0.0);
  log_debug("%s: Setting plot frequency: %g.", model->name, T);
  model->plot_every = T;
}

void model_save_every(model_t* model, int n)
{
  ASSERT(n > 0);
  log_debug("%s: Setting save frequency: %d.", model->name, n);
  model->save_every = n;
}

void model_load_from(model_t* model, int step)
{
  ASSERT(step >= 0);
  model->load_step = step;
}

void model_set_diagnostic_mode(model_t* model, model_diag_mode_t mode)
{
  model->diag_mode = mode;
}

model_vtable model_get_vtable(model_t* model)
{
  return model->vtable;
}

void model_set_vtable(model_t* model, model_vtable vtable)
{
  model->vtable = vtable;
}

