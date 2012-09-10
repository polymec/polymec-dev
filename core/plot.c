#include <stdlib.h>
#include <string.h>
#include "core/plot.h"

#ifdef __cplusplus
extern "C" {
#endif

struct plot_interface_t
{
  void* context;
  char* name;
  plot_vtable vtable;
};

plot_interface_t* plot_interface_new(void* context, const char* name, plot_vtable vtable)
{
  ASSERT(name != NULL);
  plot_interface_t* plot = malloc(sizeof(plot_interface_t));
  plot->context = context;
  plot->name = strdup(name);
  plot->vtable = vtable;
  return plot;
}

void plot_free(plot_interface_t* plot)
{
  if ((plot->context != NULL) && (plot->vtable.dtor != NULL))
    plot->vtable.dtor(plot->context);
  free(plot->name);
  free(plot);
}

char* plot_name(plot_interface_t* plot)
{
  return plot->name;
}

void plot_time_series(plot_interface_t* plot, 
                      const char* title,
                      double* data,
                      int num_data)
{
  plot->vtable.plot_time_series(plot->context, title, data, num_data);
}

void plot_mesh(plot_interface_t* plot, 
               mesh_t* mesh)
{
  plot->vtable.plot_mesh(plot->context, mesh);
}

void plot_mesh_field(plot_interface_t* plot, 
                     const char* title,
                     mesh_centering_t centering,
                     double* data)
{
  plot->vtable.plot_mesh_field(plot->context, title, centering, data);
}

void plot_lite_mesh(plot_interface_t* plot, 
                    lite_mesh_t* lite_mesh)
{
  plot->vtable.plot_lite_mesh(plot->context, lite_mesh);
}

void plot_lite_mesh_field(plot_interface_t* plot, 
                          const char* title,
                          lite_mesh_centering_t centering,
                          double* data)
{
  plot->vtable.plot_lite_mesh_field(plot->context, title, centering, data);
}

void plot_flush(plot_interface_t* plot)
{
  plot->vtable.flush(plot->context);
}

#ifdef __cplusplus
}
#endif

