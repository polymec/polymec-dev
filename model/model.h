// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_MODEL_H
#define POLYMEC_MODEL_H

#include "core/polymec.h"
#include "core/options.h"
#include "core/st_func.h"
#include "model/probe.h"

/// \addtogroup model model
///@{

/// \class model
/// A model is a numerical model of a physical phenomenon.
typedef struct model_t model_t;

/// \struct model_vtable
/// This virtual table must be implemented by any model.
typedef struct
{
  /// A function for initializing the model at time t.
  void (*init)(void* context, real_t t);

  /// A function for calculating the maximum step size.
  real_t (*max_dt)(void* context, real_t, char* reason);

  /// A function for advancing the model from t by a step of maximum size max_dt.
  /// Returns the actual size of the time step.
  real_t (*advance)(void* context, real_t max_dt, real_t t);

  /// A function for work to be performed after a run completes.
  void (*finalize)(void* context, int step, real_t t);

  /// A function for loading the model's state. All data within a model must be
  /// loaded from the file in order to ensure that the state is exactly
  /// preserved. This function is called INSTEAD OF model_init when saved data
  /// is loaded.
  /// \returns true if the load is successful, false if not.
  bool (*load)(void* context, const char* file_prefix, const char* directory, real_t* time, int step);

  /// A function for saving the model's state to the given I/O interface.
  void (*save)(void* context, const char* file_prefix, const char* directory, real_t time, int step);

  /// A function for plotting the model to the given I/O interface.
  void (*plot)(void* context, const char* file_prefix, const char* directory, real_t time, int step);

  /// A function for performing work when a probe is added.
  void (*add_probe)(void* context, void* probe_context);

  /// A destructor function for the context object (if any).
  void (*dtor)(void* context);
} model_vtable;

/// \enum model_parallelism_t
/// This type is used to identify the degree of parallelism that a given
/// model can take advantage of. We distinguish between five different "types"
/// of parallelism, based on the implementation of the underlying model:
/// 1. MODEL_SERIAL_SINGLETON: a model that is only intended to be run serially,
///    and of which only a single instance is allowed per process, likely because
///    of the use of global variables/resources.
/// 2. MODEL_SERIAL: a model that is only intended to be run serially, but of
///    which several instances can coexist in memory.
/// 3. MODEL_MPI_SINGLETON: a model that uses MPI, of which only a single
///    instance is allowed per process, likely because of the use of global
///    variables/resources.
/// 4. MODEL_MPI: a model that uses MPI, of which several instances can coexist
///    within a process, likely running different physics. Models of this type
///    are not assumed to be thread-safe, however.
/// 5. MODEL_MPI_THREAD_SAFE: a model that uses MPI and is also thread-safe,
///    so that several threads can be active within model methods simultaneously.
typedef enum
{
  MODEL_SERIAL_SINGLETON,
  MODEL_SERIAL,
  MODEL_MPI_SINGLETON,
  MODEL_MPI,
  MODEL_MPI_THREAD_SAFE
} model_parallelism_t;

/// \enum model_diag_mode_t
/// This type identifies one of two possible modes of operation for
/// diagnostics: plots and probe measurements:
/// 1. MODEL_DIAG_EXACT_TIME: This tells the model to create plots and
///    acquire probe measurements at the exact time at which they are requested,
///    altering the simulation to set time step sizes so that these data can
///    be collected.
/// 2. MODEL_DIAG_NEAREST_STEP: This tells the model to create plots and
///    acquire probe measurements on the step that falls nearest its requested
///    time. In this setting, time steps are not chosen to collect measurements.
typedef enum
{
   MODEL_DIAG_EXACT_TIME,
   MODEL_DIAG_NEAREST_STEP
} model_diag_mode_t;

/// Creates an instance of a model with the given name and characteristics.
/// The name should uniquely identify the model, BUT SHOULD NOT BE SPECIFIC
/// TO THE INSTANCE OF THE MODEL. Constructors that use this function should
/// return objects that can be safely destroyed with model_free below. This
/// means that all data members should be properly initialized in a way that
/// they can be destroyed.
/// \memberof model
model_t* model_new(const char* name,
                   void* context,
                   model_vtable vtable,
                   model_parallelism_t parallelism);

/// Destroys the model.
/// \memberof model
void model_free(model_t* model);

/// Returns the name of the model.
/// \memberof model
char* model_name(model_t* model);

/// Returns the context object associated with the model (if any).
/// \memberof model
void* model_context(model_t* model);

/// Allows iteration over all available model probe data. Set *pos to 0 to
/// reset the iteration.
/// \memberof model
bool model_next_probe_data(model_t* model,
                           int* pos,
                           char** quantity,
                           probe_data_array_t** data);

/// Returns the data array for the given probed quantity in the model, or
/// NULL if no such quantity is tracked.
/// \memberof model
probe_data_array_t* model_probe_data(model_t* model, const char* quantity);

/// Returns the degree of parallelism supported by this model.
/// \memberof model
model_parallelism_t model_parallelism(model_t* model);

/// Sets a flag to decide whether this model intercepts the SIGINT and SIGTERM
/// signals to make model steps (\ref model_advance) uninterruptible. By default,
/// a model does intercept these signals. This helps with making simulations
/// more robust, but can cause problems, say, when debugging simulations with
/// infinite loops.
/// \param [in] flag If true, the model will intercept SIGINT and SIGTERM,
///                  interrupting a simulation after the completion of a step.
///                  If false, the model won't intercept these signals.
/// \memberof model
void model_handle_signals(model_t* model, bool flag);

/// Adds a probe to the model which will acquire specific data at each of the
/// specified times. The model assumes ownership of the probe but copies the
/// array of times.
/// \memberof model
void model_add_probe(model_t* model,
                     probe_t* probe,
                     real_t* acq_times,
                     size_t num_acq_times);

/// Initializes the model at the given time.
/// \memberof model
void model_init(model_t* model, real_t t);

/// Loads the model's state. This is called instead of model_init when
/// input instructs the model to load from a given step. Returns true if
/// the model was successfully loaded, false if not.
/// \memberof model
bool model_load(model_t* model, int step);

/// Sets the initial time step to be taken by the model.
/// \memberof model
void model_set_initial_dt(model_t* model, real_t dt0);

/// Returns the initial time step to be taken by the model
/// (after initialization only).
/// \memberof model
real_t model_initial_dt(model_t* model);

/// Sets the largest permissible (positive) time step that can be taken by the
/// model.
/// \memberof model
void model_set_max_dt(model_t* model, real_t max_dt);

/// Returns the largest permissible time step that can be taken by the model.
/// \param [out] reason If non-NULL, this stores an internal string giving the
///                     reason for the selected maximum time step.
/// \memberof model
real_t model_max_dt(model_t* model, char** reason);

/// Sets the smallest permissible (non-negative) time step that can be taken by the
/// model, below which a simulation will be terminated.
/// \memberof model
void model_set_min_dt(model_t* model, real_t min_dt);

// Returns the smallest permissible time step that can be taken by the model,
// below which a simulation will be terminated.
/// \memberof model
real_t model_min_dt(model_t* model);

/// Advances the model by a single time step of maximum size max_dt.
/// If signal handling is enabled for the model (by default or via \ref model_handle_signals),
/// the model intercepts SIGINT and SIGTERM signals so that this step is
/// uninterruptible.
/// \param [in] max_dt The maximum size of the step to take.
/// \returns the size of the actual time step.
/// \memberof model
real_t model_advance(model_t* model, real_t max_dt);

/// Returns the current simulation time for the model.
/// \memberof model
real_t model_time(model_t* model);

/// Returns the current step for the model.
/// \memberof model
int model_step(model_t* model);

/// Performs any post-simulation work for the model.
/// \memberof model
void model_finalize(model_t* model);

/// Saves the model's state.
/// \memberof model
void model_save(model_t* model);

/// Plots the model's state.
/// \memberof model
void model_plot(model_t* model);

/// Acquires data from any probes added to the model, at the current time.
/// \memberof model
void model_acquire(model_t* model);

/// Runs a simulation of the model from time t1 to t2, or for a maximum of
/// max_steps. Sets up a signal handler to intercept SIGINT so that a single
/// step (a call to \ref model_advance) is uninterruptible.
/// \param [in] t1 The start time for the run.
/// \param [in] t2 The end time for the run. Must be greater than t1.
/// \param [in] max_steps The maximum number of steps in the run.
/// \memberof model
void model_run(model_t* model, real_t t1, real_t t2, int max_steps);

/// Sets the name of the prefix for naming simulation data files.
/// \memberof model
void model_set_prefix(model_t* model, const char* prefix);

/// Sets the absolute path in which the simulation will be run.
/// \memberof model
void model_set_directory(model_t* model, const char* directory);

/// Tells the model to plot a result every T time units.
/// \memberof model
void model_plot_every(model_t* model, real_t T);

/// Tells the model to plot a result every n steps.
/// \memberof model
void model_save_every(model_t* model, int n);

/// Tells the model to load its state from a given step.
/// \memberof model
void model_load_from(model_t* model, int step);

/// Sets the diagnostic mode for the model to collect measurements.
/// \memberof model
void model_set_diagnostic_mode(model_t* model, model_diag_mode_t mode);

/// Retrieves (a copy of) the virtual table for the given model.
/// \memberof model
model_vtable model_get_vtable(model_t* model);

/// Overrides the given model's virtual table with the given one. Be careful
/// when using this!
/// \memberof model
void model_set_vtable(model_t* model, model_vtable vtable);

///@}

#endif

