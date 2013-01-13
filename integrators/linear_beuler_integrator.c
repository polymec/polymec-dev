#include "integrators/linear_beuler_integrator.h"
#include "sundials/sundials_spgmr.h"
#include "sundials/sundials_spbcgs.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
  // Context stuff.
  void* context;
  integrator_dtor dtor;

  // Solver stuff.
  void* solver;
  ATimesFn Ax;
  linear_beuler_compute_rhs compute_rhs;
  int precond_type;
  PSolveFn precond;
  int gram_schmidt;
  double delta;
  int max_restarts;
  N_Vector x, b;

  // Solver metadata.
  double res_norm;
  int nli, nps;
} linear_beuler_t;

static void gmres_init(void* context, double t, double* solution, int N)
{
}

static void gmres_step(void* context, double t1, double t2, double* solution, int N)
{
}

static void gmres_dtor(void* context)
{
}

integrator_t* gmres_linear_beuler_integrator_new(void* context, 
                                                 ATimesFn Ax,
                                                 linear_beuler_compute_rhs compute_rhs,
                                                 int precond_type,
                                                 PSolveFn precond,
                                                 int gram_schmidt,
                                                 double delta,
                                                 int max_restarts,
                                                 integrator_dtor dtor)
{
  ASSERT(Ax != NULL);
  ASSERT(compute_rhs != NULL);
  ASSERT((precond_type == PREC_NONE) || (precond_type == PREC_LEFT) ||
         (precond_type == PREC_RIGHT) || (precond_type == PREC_BOTH));
  ASSERT((precond != NULL) || (precond_type == PREC_NONE));
  ASSERT((gram_schmidt == CLASSICAL_GS) || (gram_schmidt == MODIFIED_GS));
  ASSERT(delta > 0.0);
  ASSERT(max_restarts >= 0);
  linear_beuler_t* beuler = malloc(sizeof(linear_beuler_t));
  beuler->context = context;
  beuler->solver = NULL;
  beuler->Ax = Ax;
  beuler->compute_rhs = compute_rhs;
  beuler->precond_type = precond_type;
  beuler->precond = precond;
  beuler->gram_schmidt = gram_schmidt;
  beuler->delta = delta;
  beuler->max_restarts = max_restarts;
  beuler->dtor = dtor;
  beuler->res_norm = 0.0;
  beuler->nli = 0;
  beuler->nps = 0;
  integrator_vtable vtable = { .init = gmres_init, .step = gmres_step, .dtor = gmres_dtor };
  return integrator_new("Backward Euler (GMRES)", beuler, vtable, 1, INTEGRATOR_IMPLICIT);
}

static void bicgstab_init(void* context, double t, double* solution, int N)
{
}

static void bicgstab_step(void* context, double t1, double t2, double* solution, int N)
{
}

static void bicgstab_dtor(void* context)
{
}

integrator_t* bicgstab_linear_beuler_integrator_new(void* context, 
                                                    ATimesFn Ax,
                                                    linear_beuler_compute_rhs compute_rhs,
                                                    int precond_type,
                                                    PSolveFn precond,
                                                    double delta,
                                                    integrator_dtor dtor)
{
  ASSERT(Ax != NULL);
  ASSERT(compute_rhs != NULL);
  ASSERT((precond_type == PREC_NONE) || (precond_type == PREC_LEFT) ||
         (precond_type == PREC_RIGHT) || (precond_type == PREC_BOTH));
  ASSERT((precond != NULL) || (precond_type == PREC_NONE));
  ASSERT(delta > 0.0);
  linear_beuler_t* beuler = malloc(sizeof(linear_beuler_t));
  beuler->context = context;
  beuler->solver = NULL;
  beuler->Ax = Ax;
  beuler->compute_rhs = compute_rhs;
  beuler->precond_type = precond_type;
  beuler->precond = precond;
  beuler->delta = delta;
  beuler->dtor = dtor;
  beuler->res_norm = 0.0;
  beuler->nli = 0;
  beuler->nps = 0;
  integrator_vtable vtable = { .init = bicgstab_init, .step = bicgstab_step, .dtor = bicgstab_dtor };
  return integrator_new("Backward Euler (GMRES)", beuler, vtable, 1, INTEGRATOR_IMPLICIT);
}

void linear_beuler_integrator_get_info(integrator_t* integrator,
                                       double* res_l2_norm,
                                       int* num_linear_iterations,
                                       int* num_precond_solves)
{
  linear_beuler_t* beuler = integrator_context(integrator);
  *res_l2_norm = beuler->res_norm;
  *num_linear_iterations = beuler->nli;
  *num_precond_solves = beuler->nps;
}

#ifdef __cplusplus
}
#endif

