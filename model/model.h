// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_MODEL_H
#define POLYMEC_MODEL_H

#include "core/polymec.h"
#include "core/options.h"
#include "model/interpreter.h"

// The maximum amount of storage allowed for an explanation of the 
// time step choice.
#define POLYMEC_MODEL_MAXDT_REASON_SIZE 2048

// A model is a numerical model of a physical phenomenon.
typedef struct model_t model_t;

// A model constructor function for creating an object context.
typedef model_t* (*model_ctor)();

// A function for setting the MPI communicator to be used by the model. This 
// communicator will serve the role of MPI_COMM_WORLD in parallel communications 
// within the model.
typedef void (*model_set_comm_func)(void* context, MPI_Comm comm);

// A function for reading input from an interpreter into the model.
typedef void (*model_read_input_func)(void* context, interpreter_t* interpreter, options_t* options);

// A function for reading custom (non-Lua) input from a file into the model.
typedef void (*model_read_custom_input_func)(void* context, const char* input_string, options_t* options);

// A function for initializing the model at time t.
typedef void (*model_init_func)(void* context, real_t t);

// A function for calculating the maximum step size.
typedef real_t (*model_max_dt_func)(void* context, real_t t, char* reason);

// A function for advancing the model from t by a step of maximum size max_dt.
// Returns the actual size of the time step.
typedef real_t (*model_advance_func)(void* context, real_t max_dt, real_t t);

// A function for work to be performed after a run completes.
typedef void (*model_finalize_func)(void* context, int step, real_t t);

// A function for loading the model's state. All data within a model must be 
// loaded from the file in order to ensure that the state is exactly 
// preserved. This function is called INSTEAD OF model_init when saved data
// is loaded.
typedef void (*model_load_func)(void* context, const char* file_prefix, const char* directory, real_t* time, int step);

// A function for saving the model's state to the given I/O interface.
typedef void (*model_save_func)(void* context, const char* file_prefix, const char* directory, real_t time, int step);

// A function for plotting the model to the given I/O interface.
typedef void (*model_plot_func)(void* context, const char* file_prefix, const char* directory, real_t time, int step);

// A function for computing error norms for the computed solution, compared 
// with an analytic solution. The error norms can be specific to the model.
typedef void (*model_compute_error_norms_func)(void* context, st_func_t* solution, real_t t, real_t* norms);

// A destructor function for the context object (if any).
typedef void (*model_dtor)(void* context);

// This virtual table must be implemented by any model.
typedef struct 
{
  model_set_comm_func            set_comm;
  model_read_input_func          read_input;
  model_read_custom_input_func   read_custom_input;
  model_init_func                init;
  model_max_dt_func              max_dt;
  model_advance_func             advance;
  model_finalize_func            finalize;
  model_load_func                load;
  model_save_func                save;
  model_plot_func                plot;
  model_compute_error_norms_func compute_error_norms;
  model_dtor                     dtor;
} model_vtable;

// This is used to construct a table of models for use with 
// multi_model_main(), below.
typedef struct 
{
  char* model_name;
  model_ctor model_constructor;
} model_dispatch_t;

// This is used to terminate the table passed to multi_model_main(), below.
static const model_dispatch_t END_OF_MODELS = {(char*)"END", NULL};

// This type is used to identify the degree of parallelism that a given 
// model can take advantage of. We distinguish between five different "types"
// of parallelism, based on the implementation of the underlying model:
// 1. MODEL_SERIAL_SINGLETON: a model that is only intended to be run serially,
//    and of which only a single instance is allowed per process, likely because
//    of the use of global variables/resources.
// 2. MODEL_SERIAL: a model that is only intended to be run serially, but of 
//    which several instances can coexist in memory.
// 3. MODEL_MPI_SINGLETON: a model that uses MPI, of which only a single 
//    instance is allowed per process, likely because of the use of global 
//    variables/resources.
// 4. MODEL_MPI: a model that uses MPI, of which several instances can coexist
//    within a process, likely running different physics. Models of this type 
//    are not assumed to be thread-safe, however.
// 5. MODEL_MPI_THREAD_SAFE: a model that uses MPI and is also thread-safe, 
//    so that several threads can be active within model methods simultaneously.
typedef enum
{
  MODEL_SERIAL_SINGLETON, 
  MODEL_SERIAL, 
  MODEL_MPI_SINGLETON,
  MODEL_MPI, 
  MODEL_MPI_THREAD_SAFE 
} model_parallelism_t;

// Creates an instance of a model with the given name and characteristics.
// The name should uniquely identify the model, BUT SHOULD NOT BE SPECIFIC 
// TO THE INSTANCE OF THE MODEL. Constructors that use this function should 
// return objects that can be safely destroyed with model_free below. This 
// means that all data members should be properly initialized in a way that 
// they can be destroyed.
model_t* model_new(const char* name, 
                   void* context, 
                   model_vtable vtable, 
                   docstring_t* doc,
                   model_parallelism_t parallelism);

// Destroys the model.
void model_free(model_t* model);

// Returns the name of the model.
char* model_name(model_t* model);

// Returns the context object associated with the model (if any).
void* model_context(model_t* model);

// Returns an internal pointer to the interpreter that the model uses to 
// parse input files.
interpreter_t* model_interpreter(model_t* model);

// Enables an interpreter for the model with the given set of variables
// to validate against types.
void model_enable_interpreter(model_t* model, interpreter_validation_t* valid_inputs);

// Returns the degree of parallelism supported by this model.
model_parallelism_t model_parallelism(model_t* model);

// Prints usage information for the model to the given file stream.
void model_usage(model_t* model, FILE* stream);

// This is the type of function that needs to be implemented to 
// run a benchmark calculation.
typedef void (*model_benchmark_function_t)();

// Registers the given benchmark name, function, and description with 
// this model. The description is consumed if given.
void model_register_benchmark(model_t* model, const char* benchmark, model_benchmark_function_t function, docstring_t* description);

// Writes a description of the given benchmark to the given file stream.
void model_describe_benchmark(model_t* model, const char* benchmark, FILE* stream);

// Runs the given benchmark problem for the model.
void model_run_benchmark(model_t* model, const char* benchmark);

// Runs all benchmark problems for the model.
void model_run_all_benchmarks(model_t* model);

// Tells the model to read the contents of the given string as input, 
// with values overridden by options.
void model_read_input_string(model_t* model, const char* input);

// Tells the model to read the contents of the file with the given name 
// as input, with values overwritten by options.
void model_read_input_file(model_t* model, const char* file);

// Defines a point observation by its name, number of components, and 
// calculation function, and at the given point.
void model_define_point_observation(model_t* model, 
                                    const char* name, 
                                    real_t (*compute_point_observation)(void* context, 
                                                                        point_t* x,
                                                                        real_t t),
                                    point_t* point);

// Defines a global (non-point-specific) observation by its name, 
// number of components, and calculation function. Examples of global 
// observations are integrated quantities, maxima, and minima.
void model_define_global_observation(model_t* model, 
                                     const char* name, 
                                     real_t (*compute_global_observation)(void* context, 
                                                                          real_t t));

// Adds the given observation to the set that will be recorded by the 
// model during a simulation.
void model_observe(model_t* model, const char* observation);

// Sets the observation times for the model, overwriting any old observation
// times.
void model_set_observation_times(model_t* model, real_t* times, int num_times);

// Initializes the model at the given time.
void model_init(model_t* model, real_t t);

// Loads the model's state. This is called instead of model_init when 
// input instructs the model to load from a given step.
void model_load(model_t* model, int step);

// Sets the initial time step to be taken by the model.
void model_set_initial_dt(model_t* model, real_t dt0);

// Returns the initial time step to be taken by the model 
// (after initialization only).
real_t model_initial_dt(model_t* model);

// Sets the largest permissible (positive) time step that can be taken by the 
// model.
void model_set_max_dt(model_t* model, real_t max_dt);

// Returns the largest permissible time step that can be taken by the model.
real_t model_max_dt(model_t* model, char* reason);

// Advances the model by a single time step of maximum size max_dt, returning 
// the size of the actual time step.
real_t model_advance(model_t* model, real_t max_dt);

// Returns the current simulation time for the model.
real_t model_time(model_t* model);

// Returns the current step for the model.
int model_step(model_t* model);

// Performs any post-simulation work for the model.
void model_finalize(model_t* model);

// Saves the model's state.
void model_save(model_t* model);

// Plots the model's state.
void model_plot(model_t* model);

// Records any desired observations for the model at the current time.
void model_record_observations(model_t* model);

// Given a model with a computed solution, compute the error norms for 
// the solution versus a specified analytic solution (at the given time).
void model_compute_error_norms(model_t* model, st_func_t* solution, real_t* error_norms);

// Runs a simulation of the model from time t1 to t2, or for a maximum of 
// max_steps.
void model_run(model_t* model, real_t t1, real_t t2, int max_steps);

// Runs a batch of simulations whose inputs are found in files whose names 
// appear in input_files. The simulations will be run concurrently to the 
// degree possible. The input files will only be read by process 0 of the 
// MPI_COMM_WORLD communicator.
void model_run_files(model_t* model, 
                     char** input_files,
                     size_t num_input_files);

// Sets the name of the simulation within the model. This name 
// will be used to identify and/or generate names of plot and save files.
void model_set_sim_name(model_t* model, const char* sim_name);

// Sets the absolute path in which the simulation will be run. All output
// will be generated relative to this path.
void model_set_sim_path(model_t* model, const char* sim_path);

// This function implements a simple driver for a model and behaves
// in the same way as a main() function, returning 0 on success and nonzero
// on failure. Arguments:
// model_name  - The name of the model as it is called by the user.
// constructor - A constructor function that uses a given set of option 
//               strings to construct a model object.
// argc        - The number of command line arguments.
// argv        - The command line arguments.
int model_main(const char* model_name, model_ctor constructor, int argc, char* argv[]);

// This function implements a driver for a collection of models and behaves
// in the same way as model_main(), except that all of its commands take an 
// extra argument that identifies one of several models from the given 
// model dispatch table.
// Arguments:
// model_table - An array of key-value pairs mapping names of models (strings)
//               to model constructors, and terminated with END_OF_MODELS.
// argc        - The number of command line arguments.
// argv        - The command line arguments.
int multi_model_main(model_dispatch_t model_table[], 
                     int argc, 
                     char* argv[]);

// This version of model_main() does not support "commands," and only runs 
// a simulation given an input file in the form:
// <executable> <input> [options] ...
int model_minimal_main(const char* model_name, model_ctor constructor, int argc, char* argv[]);

// Use this to report a convergence rate from within a benchmark. This can be 
// used to determine whether the given benchmark "passed."
void model_report_conv_rate(real_t conv_rate, real_t sigma);

// Use this to report an error norm from within a benchmark. This can be 
// used to determine whether the given benchmark "passed."
void model_report_error_norm(real_t error_norm);

#endif

