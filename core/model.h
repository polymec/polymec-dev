#ifndef ARBI_MODEL_H
#define ARBI_MODEL_H

#include "core/arbi.h"
#include "core/options.h"
#include "core/io.h"
#include "core/plot.h"

#ifdef __cplusplus
extern "C" {
#endif

// The maximum amount of storage allowed for an explanation of the 
// time step choice.
#define ARBI_MODEL_MAXDT_REASON_SIZE 2048

// A model is a numerical model of a physical phenomenon.
typedef struct model_t model_t;

// A model constructor function for creating an object context.
typedef void* (*model_ctor)(options_t* options);

// A function for running a benchmark calculation.
typedef void (*model_run_bench_func)(const char*);

// A function for initializing the model.
typedef void (*model_init_func)(void*, double);

// A function for calculating the maximum step size.
typedef double (*model_max_dt_func)(void*, double, char*);

// A function for advancing the model.
typedef void (*model_advance_func)(void*, double, double);

// A function for loading the model's state from an I/O interface.
typedef void (*model_load_func)(void*, io_interface_t*, double*, int*);

// A function for dumping the model's state to an I/O interface.
typedef void (*model_dump_func)(void*, io_interface_t*, double, int);

// A function for plotting the model to a plot interface.
typedef void (*model_plot_func)(void*, plot_interface_t*, double, int);

// A destructor function for the context object (if any).
typedef void (*model_dtor)(void*);

// This virtual table must be implemented by any model.
typedef struct 
{
  model_run_bench_func  run_benchmark;
  model_init_func       init;
  model_max_dt_func     max_dt;
  model_advance_func    advance;
  model_load_func       load;
  model_dump_func       dump;
  model_plot_func       plot;
  model_dtor            dtor;
} model_vtable;

// Creates an instance a model with the given name and characteristics.
model_t* model_new(const char* name, void* context, model_vtable vtable);

// Destroy the model.
void model_free(model_t* model);

// Returns the name of the model.
char* model_name(model_t* model);

// Returns the context object associated with the model (if any).
void* model_context(model_t* model);

// Print usage information for the model to the given file stream.
void model_usage(model_t* model, FILE* stream);

// Associate the given list of benchmarks with this model.
// Here, benchmarks is a NULL-terminated array of strings.
void model_register_benchmarks(model_t* model, const char** benchmarks);

// Returns an array containing the names of benchmarks, plus its length.
void model_get_benchmarks(model_t* model, char*** benchmarks, int* num_benchmarks);

// Runs the given benchmark problem for the model.
void model_run_benchmark(model_t* model, const char* benchmark);

// Runs all benchmark problems for the model.
void model_run_benchmarks(model_t* model);

// Initialize the model at the given time.
void model_init(model_t* model, double t);

// Returns the largest permissible time step that can be taken by the model
// starting at time t.
double model_max_dt(model_t* model, double t, char* reason);

// Advance the model by a single time step of size dt, beginning at time t.
void model_advance(model_t* model, double t, double dt);

// Load the model's state from the given I/O interface.
void model_load(model_t* model, io_interface_t* io, double* t, int* step);

// Dump the model's state to the given I/O interface.
void model_dump(model_t* model, io_interface_t* io, double t, int step);

// Plot the model's state to the given plot interface.
void model_plot(model_t* model, plot_interface_t* plot, double t, int step);

// Run a simulation of the model from time t1 to t2.
void model_run(model_t* model, double t1, double t2);

#ifdef __cplusplus
}
#endif

#endif

