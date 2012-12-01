#ifndef POLYMEC_MODEL_H
#define POLYMEC_MODEL_H

#include "core/polymec.h"
#include "core/options.h"
#include "core/interpreter.h"
#include "core/io.h"

#ifdef __cplusplus
extern "C" {
#endif

// The maximum amount of storage allowed for an explanation of the 
// time step choice.
#define POLYMEC_MODEL_MAXDT_REASON_SIZE 2048

// A model is a numerical model of a physical phenomenon.
typedef struct model_t model_t;

// A model constructor function for creating an object context.
typedef model_t* (*model_ctor)(options_t*);

// A function for reading input from an interpreter into the model.
typedef void (*model_read_input_func)(void*, interpreter_t*, options_t*);

// A function for initializing the model.
typedef void (*model_init_func)(void*, double);

// A function for calculating the maximum step size.
typedef double (*model_max_dt_func)(void*, double, char*);

// A function for advancing the model.
typedef void (*model_advance_func)(void*, double, double);

// A function for work to be performed after a run completes.
typedef void (*model_finalize_func)(void*, int, double);

// A function for loading the model's state from the given I/O interface.
typedef void (*model_load_func)(void*, io_interface_t*, double*, int);

// A function for saving the model's state to the given I/O interface.
typedef void (*model_save_func)(void*, io_interface_t*, double, int);

// A function for plotting the model to the given I/O interface.
typedef void (*model_plot_func)(void*, io_interface_t*, double, int);

// A function for computing error norms for the computed solution, compared 
// with an analytic solution. The error norms can be specific to the model.
typedef void (*model_compute_error_norms_func)(void*, st_func_t*, double, double*);

// A destructor function for the context object (if any).
typedef void (*model_dtor)(void*);

// This virtual table must be implemented by any model.
typedef struct 
{
  model_read_input_func          read_input;
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

// Creates an instance a model with the given name and characteristics, 
// along with any relevant options.
model_t* model_new(const char* name, void* context, model_vtable vtable, options_t* options);

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
typedef void (*model_benchmark_function_t)(options_t*);

// Registers the given benchmark name, function, and description with 
// this model.
void model_register_benchmark(model_t* model, const char* benchmark, model_benchmark_function_t function, const char* description);

// Runs the given benchmark problem for the model.
void model_run_benchmark(model_t* model, const char* benchmark, options_t* options);

// Runs all benchmark problems for the model.
void model_run_all_benchmarks(model_t* model, options_t* options);

// Tells the model to read the contents of the given string as input, 
// tempered by options.
void model_read_input_string(model_t* model, const char* input, options_t* options);

// Tells the model to read the contents of the given file as input,
// tempered by options.
void model_read_input_file(model_t* model, const char* file, options_t* options);

// Initializes the model at the given time.
void model_init(model_t* model, double t);

// Returns the largest permissible time step that can be taken by the model
// starting at time t.
double model_max_dt(model_t* model, char* reason);

// Advances the model by a single time step of size dt.
void model_advance(model_t* model, double dt);

// Performs any post-simulation work for the model.
void model_finalize(model_t* model);

// Loads the model's state from its I/O interface.
void model_load(model_t* model, int step);

// Saves the model's state to its I/O interface.
void model_save(model_t* model);

// Plots the model's state to its I/O interface.
void model_plot(model_t* model);

// Given a model with a computed solution, compute the error norms for 
// the solution versus a specified analytic solution (at the given time).
void model_compute_error_norms(model_t* model, st_func_t* solution, double* error_norms);

// Runs a simulation of the model from time t1 to t2, or for a maximum of 
// max_steps.
void model_run(model_t* model, double t1, double t2, int max_steps);

// Sets the name of the simulation within the model. This name 
// will be used to identify and/or generate names of plot and save files.
void model_set_sim_name(model_t* model, const char* sim_name);

// Sets the I/O interface that will be used to save the model's state.
void model_set_saver(model_t* model, io_interface_t* saver);

// Sets the I/O interface that will be used to plot the model's data.
void model_set_plotter(model_t* model, io_interface_t* plotter);

// This function implements a simple driver for a model and behaves
// in the same way as a main() function, returning 0 on success and nonzero
// on failure. Arguments:
// model_name  - The name of the model as it is called by the user.
// constructor - A constructor function that uses a given set of option 
//               strings to construct a model object.
// argc        - The number of command line arguments.
// argv        - The command line arguments.
int model_main(const char* model_name, model_ctor constructor, int argc, char* argv[]);

#ifdef __cplusplus
}
#endif

#endif

