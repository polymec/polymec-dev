#include "integrators/ark_integrator.h"

typedef struct
{
  integrator_t *ex_int, *imp_int;
} ark_t;

static void asirk1a_init(void* context, int N)
{
  ark_t* ark = context;
  integrator_init(ark->ex_int, N);
  integrator_init(ark->imp_int, N);
}

static void asirk1a_step(void* context, double t1, double t2, double* solution, int N)
{
  static const double omega1 = 1.0, a1 = 1.0;
  ASSERT(t2 > t1);
  double h = t2 - t1;
  ark_t* ark = context;
}

static void asirk1b_init(void* context, int N)
{
  ark_t* ark = context;
  integrator_init(ark->ex_int, N);
  integrator_init(ark->imp_int, N);
}

static void asirk1b_step(void* context, double t1, double t2, double* solution, int N)
{
  ASSERT(t2 > t1);
  double h = t2 - t1;
  ark_t* ark = context;
}

static void asirk1c_init(void* context, int N)
{
  ark_t* ark = context;
  integrator_init(ark->ex_int, N);
  integrator_init(ark->imp_int, N);
}

static void asirk1c_step(void* context, double t1, double t2, double* solution, int N)
{
  ASSERT(t2 > t1);
  double h = t2 - t1;
  ark_t* ark = context;
}

static void ark_dtor(void* context)
{
  ark_t* ark = context;
  integrator_free(ark->ex_int);
  integrator_free(ark->imp_int);
  free(ark);
}

integrator_t* ark_integrator_new(ark_integrator_type_t type,
                                 integrator_t* explicit_integrator,
                                 integrator_t* implicit_integrator)
{
  ASSERT(explicit_integrator != NULL);
  ASSERT(implicit_integrator != NULL);

  ark_t* ark = malloc(sizeof(ark_t));
  ark->ex_int = explicit_integrator;
  ark->imp_int = implicit_integrator;
  integrator_vtable vtable;
  char name[1024];
  int order;
  switch (type)
  {
    case ASIRK_1A: order = 1;
                   vtable.init = asirk1a_init;
                   vtable.step = asirk1a_step;
                   strcpy(name, "ASIRK-1A");
                   break;
    case ASIRK_1B: order = 1;
                   vtable.init = asirk1b_init;
                   vtable.step = asirk1b_step;
                   strcpy(name, "ASIRK-1B");
                   break;
    case ASIRK_1C: order = 1; 
                   vtable.init = asirk1c_init;
                   vtable.step = asirk1c_step;
                   strcpy(name, "ASIRK-1C");
                   break;
    default: break;
  }
  vtable.dtor = ark_dtor;
  return integrator_new(name, ark, vtable, order, INTEGRATOR_SEMI_IMPLICIT);
}

