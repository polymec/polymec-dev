#include <stdlib.h>
#include <string.h>
#include <float.h>
#include "core/model.h"

#ifdef __cplusplus
extern "C" {
#endif

struct model_t 
{
  void* context;
  char* name;
  char* usage;
  model_vtable vtable;
};

model_t* model_new(void* context, 
                   const char* name, 
                   const char* usage_string,
                   model_vtable vtable,
                   options_t* options)
{
  ASSERT(name != NULL);
  ASSERT(usage_string != NULL);
  model_t* model = malloc(sizeof(model_t));
  model->context = context;
  model->name = strdup(name);
  model->usage = strdup(usage_string);
  model->vtable = vtable;
  return model;
}

void model_free(model_t* model)
{
  if ((model->context != NULL) && (model->vtable.dtor != NULL))
    model->vtable.dtor(model->context);
  free(model->name);
  free(model->usage);
  free(model);
}

char* model_name(model_t* model)
{
  return model->name;
}

// Print usage information for the model to the given file stream.
void model_usage(model_t* model, FILE* stream)
{
  fprintf(stream, "%s", model->usage);
}

void model_run_benchmark(model_t* model, char* benchmark)
{
  if (model->vtable.run_benchmark != NULL)
    model->vtable.run_benchmark(model->context, benchmark);
  else
  {
    char err[1024];
    snprintf(err, 1024, "No benchmarks are defined for the model '%s'.", model->name);
    arbi_error(err);
  }
}

// Initialize the model at the given time.
void model_init(model_t* model, double t)
{
  model->vtable.init(model->context, t);
}

// Returns the largest permissible time step that can be taken by the model
// starting at time t.
double model_max_dt(model_t* model, double t, char* reason)
{
  if (model->vtable.max_dt != NULL)
    return model->vtable.max_dt(model->context, t, reason);
  else
  {
    strcpy(reason, "No time step constraints.");
    return FLT_MAX;
  }
}

void model_advance(model_t* model, double t, double dt)
{
  model->vtable.advance(model->context, t, dt);
}

void model_load(model_t* model, io_interface_t* io, double* t, int* step)
{
  model->vtable.load(model->context, io, t, step);
}

void model_dump(model_t* model, io_interface_t* io, double t, int step)
{
  model->vtable.dump(model->context, io, t, step);
}

void model_plot(model_t* model, plot_interface_t* plot, double t, int step)
{
  model->vtable.plot(model->context, plot, t, step);
}

#ifdef __cplusplus
}
#endif

