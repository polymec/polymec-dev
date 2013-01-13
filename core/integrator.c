#include "core/integrator.h"

#ifdef __cplusplus
extern "C" {
#endif

struct integrator_t 
{
  char* name;
  void* context;
  integrator_vtable vtable;
  int order;
  integrator_type_t type;
  int dim; // Number of unknowns.
};

integrator_t* integrator_new(const char* name, 
                             void* context,
                             integrator_vtable vtable,
                             int order,
                             integrator_type_t type)
{
  ASSERT(vtable.step != NULL);
  ASSERT(order > 0);
  integrator_t* integ = malloc(sizeof(integrator_t));
  integ->name = strdup(name);
  integ->context = context;
  integ->vtable = vtable;
  integ->order = order;
  integ->type = type;
  integ->dim = -1;
  return integ;
}

void integrator_free(integrator_t* integrator)
{
  if ((integrator->context != NULL) && (integrator->vtable.dtor != NULL))
    integrator->vtable.dtor(integrator->context);
  free(integrator->name);
  free(integrator);
}

char* integrator_name(integrator_t* integrator)
{
  return integrator->name;
}

void* integrator_context(integrator_t* integrator)
{
  return integrator->context;
}

int integrator_order(integrator_t* integrator)
{
  return integrator->order;
}

integrator_type_t integrator_type(integrator_t* integrator)
{
  return integrator->type;
}

void integrator_init(integrator_t* integrator, double t, double* solution, int N)
{
  integrator->vtable.init(integrator->context, t, solution, N);
  integrator->dim = N;
}

void integrator_step(integrator_t* integrator, double t1, double t2, 
                     double* solution)
{
  integrator->vtable.step(integrator->context, t1, t2, solution, integrator->dim);
}

#ifdef __cplusplus
}
#endif

