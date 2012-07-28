#ifndef ARBI_PLOT_H
#define ARBI_PLOT_H

#include "core/arbi.h"
#include "core/mesh.h"
#include "core/lite_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

// The "plot" interface type is an opaque type for descriptors used to 
// plot data in simulations.
typedef struct plot_interface_t plot_interface_t;

// A function pointer type for plotting time series data.
typedef int (*plot_time_series_func)(void*, const char*, double*, int);

// A function pointer type for plotting field data on a mesh.
typedef int (*plot_mesh_field_func)(void*, const char*, mesh_centering_t, double*, mesh_t*);

// A function pointer type for plotting field data on a lite mesh.
typedef int (*plot_lite_mesh_field_func)(void*, const char*, lite_mesh_centering_t, double*, mesh_t*);

// A destructor function for the context object (if any).
typedef void (*plot_dtor)(void*);

// This virtual table must be implemented by any I/O interface used 
// within arbi.
typedef struct 
{
  plot_time_series_func           plot_time_series;
  plot_mesh_field_func            plot_mesh_field;
  plot_lite_mesh_field_func       plot_lite_mesh_field;
  plot_dtor                       dtor;
} plot_vtable;

// Construct a plot interface (subclass) object from the given name and 
// vtable.
plot_interface_t* plot_interface_new(void* context, const char* name, plot_vtable vtable);

// Frees the given I/O interface.
void plot_free(plot_interface_t* interface);

// Plots time series data.
void plot_time_series(plot_interface_t* interface, 
                      const char* title,
                      double* data,
                      int num_data);

// Plots a mesh field.
void plot_mesh_field(plot_interface_t* interface, 
                     const char* title,
                     mesh_centering_t centering,
                     double* data,
                     mesh_t* mesh);

// Plots a lite mesh field.
void plot_lite_mesh_field(plot_interface_t* interface, 
                          const char* title,
                          lite_mesh_centering_t centering,
                          double* data,
                          lite_mesh_t* lite_mesh);

#ifdef __cplusplus
}
#endif

#endif

