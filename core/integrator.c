#include "core/integrator.h"

struct integrator_t 
{
  char* name;
  void* context;
  integrator_vtable vtable;
  int order;
  integrator_type_t type;
  int dim; // Number of unknowns.

  N_Vector temp; // Work vector for Jacobians.
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
  integ->temp = NULL;
  return integ;
}

void integrator_free(integrator_t* integrator)
{
  if ((integrator->context != NULL) && (integrator->vtable.dtor != NULL))
    integrator->vtable.dtor(integrator->context);
  if (integrator->temp != NULL)
    N_VDestroy(integrator->temp);
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

void integrator_init(integrator_t* integrator, int N)
{
  if (integrator->temp != NULL)
  {
    N_VDestroy(integrator->temp);
    integrator->temp = NULL;
  }
  integrator->vtable.init(integrator->context, N);
  integrator->dim = N;
}

void integrator_step(integrator_t* integrator, double t1, double t2, 
                     double* solution)
{
  ASSERT(t2 > t1);
  integrator->vtable.step(integrator->context, t1, t2, solution, integrator->dim);
}

void integrator_compute_Jv(integrator_t* integrator, double t, N_Vector x, N_Vector F, N_Vector v, N_Vector Jv)
{
  ASSERT(integrator->vtable.compute_Jv != NULL);
  if (integrator->temp == NULL)
    integrator->temp = N_VClone(x);
  integrator->vtable.compute_Jv(v, Jv, t, x, F, integrator->context, integrator->temp);
}

