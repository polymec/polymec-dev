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

#ifndef POLYMEC_MODEL_H
#define POLYMEC_MODEL_H

#include "core/polymec.h"
#include "core/options.h"
#include "core/interpreter.h"

// The maximum amount of storage allowed for an explanation of the 
// time step choice.
#define POLYMEC_MODEL_MAXDT_REASON_SIZE 2048

// A model is a numerical model of a physical phenomenon.
typedef struct model_t model_t;

// A model constructor function for creating an object context.
typedef model_t* (*model_ctor)();

// A function for reading input from an interpreter into the model.
typedef void (*model_read_input_func)(void* context, interpreter_t* interpreter, options_t* options);

// A function for reading custom (non-Lua) input from a file into the model.
typedef void (*model_read_custom_input_func)(void* context, const char* input_string, options_t* options);

// A function for initializing the model.
typedef void (*model_init_func)(void* context, real_t t);

// A function for calculating the maximum step size.
typedef real_t (*model_max_dt_func)(void* context, real_t t, char* reason);

// A function for advancing the model from t by a step of maximum size max_dt.
// Returns the actual size of the time step.
typedef real_t (*model_advance_func)(void* context, real_t max_dt, real_t t);

// A function for work to be performed after a run completes.
typedef void (*model_finalize_func)(void* context, int step, real_t t);

// A function for loading the model's state.
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

// Creates an instance of a model with the given name and characteristics.
model_t* model_new(const char* name, void* context, model_vtable vtable, docstring_t* doc);

// Destroys the model.
void model_free(model_t* model);

// Returns the name of the model.
char* model_name(model_t* model);

// Returns the context object associated with the model (if any).
void* model_context(model_t* model);

// The interpreter that the model uses to parse input files.
interpreter_t* model_interpreter(model_t* model);

// Enables an interpreter for the model with the given set of variables
// to validate against types.
void model_enable_interpreter(model_t* model, interpreter_validation_t* valid_inputs);

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

// Returns the largest permissible time step that can be taken by the model
// starting at time t.
real_t model_max_dt(model_t* model, char* reason);

// Advances the model by a single time step of maximum size max_dt, returning 
// the size of the actual time step.
real_t model_advance(model_t* model, real_t max_dt);

// Performs any post-simulation work for the model.
void model_finalize(model_t* model);

// Loads the model's state.
void model_load(model_t* model, int step);

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

// Sets the name of the simulation within the model. This name 
// will be used to identify and/or generate names of plot and save files.
void model_set_sim_name(model_t* model, const char* sim_name);

// This function implements a simple driver for a model and behaves
// in the same way as a main() function, returning 0 on success and nonzero
// on failure. Arguments:
// model_name  - The name of the model as it is called by the user.
// constructor - A constructor function that uses a given set of option 
//               strings to construct a model object.
// argc        - The number of command line arguments.
// argv        - The command line arguments.
int model_main(const char* model_name, model_ctor constructor, int argc, char* argv[]);

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

