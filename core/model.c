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
  model_vtable vtable;
  char** benchmarks;
  int num_benchmarks;
};

model_t* model_new(const char* name, void* context, model_vtable vtable)
{
  model_t* model = malloc(sizeof(model_t));
  model->vtable = vtable;
  model->context = context;
  model->name = strdup(name);
  model->num_benchmarks = 0;
  model->benchmarks = NULL;
  return model;
}

void model_free(model_t* model)
{
  if ((model->context != NULL) && (model->vtable.dtor != NULL))
    model->vtable.dtor(model->context);
  free(model->name);
  model_register_benchmarks(model, NULL);
  free(model);
}

char* model_name(model_t* model)
{
  return model->name;
}

void model_register_benchmarks(model_t* model, const char** benchmarks)
{
  if (model->benchmarks != NULL)
  {
    for (int i = 0; i < model->num_benchmarks; ++i)
      free(model->benchmarks[i]);
    free(model->benchmarks);
  }
  model->num_benchmarks = 0;
  char* p = (char*)benchmarks[0];
  while (p != NULL) 
    model->num_benchmarks++;
  model->benchmarks = malloc(sizeof(char*)*model->num_benchmarks);
  for (int i = 0; i < model->num_benchmarks; ++i)
    model->benchmarks[i] = strdup(benchmarks[i]);
}

void model_get_benchmarks(model_t* model, char*** benchmarks, int* num_benchmarks)
{
  *benchmarks = model->benchmarks;
  *num_benchmarks = model->num_benchmarks;
}

void model_run_benchmarks(model_t* model)
{
  for (int i = 0; i < model->num_benchmarks; ++i)
    model_run_benchmark(model, (const char*)model->benchmarks[i]);
}

void model_run_benchmark(model_t* model, const char* benchmark)
{
  if (model->vtable.run_benchmark != NULL)
    model->vtable.run_benchmark(benchmark);
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

void model_run(model_t* model, double t1, double t2)
{
  double t = t1;
  model_init(model, t);
  while (t < t2)
  {
    char reason[ARBI_MODEL_MAXDT_REASON_SIZE];
    double dt = model_max_dt(model, t, reason);
    model_advance(model, t, dt);
    t += dt;
  }
}

void* model_context(model_t* model)
{
  return model->context;
}

#ifdef __cplusplus
}
#endif

