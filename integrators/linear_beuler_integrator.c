#include "integrators/linear_beuler_integrator.h"
#include "sundials/sundials_spgmr.h"
#include "sundials/sundials_spbcgs.h"
#include "sundials/sundials_nvector.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
  // Context stuff.
  void* context;
  integrator_dtor dtor;

  // Solver stuff.
  MPI_Comm comm;
  void* solver;
  integrator_Ax_func Ax;
  integrator_compute_rhs_func compute_rhs;
  int precond_type;
  integrator_precond_func precond;
  int max_kdim;
  int gram_schmidt;
  double delta;
  int max_restarts;
  N_Vector x, b, sx, sb, s1, s2;

  // Solver metadata.
  double res_norm;
  int nli, nps;
} linear_beuler_t;

static void linear_beuler_reset(linear_beuler_t* beuler)
{
  if (beuler->solver != NULL)
  {
    SpbcgMem solver = beuler->solver;
    SpbcgFree(solver);
    N_VDestroy(beuler->x);
    N_VDestroy(beuler->b);
  }
  beuler->x = NULL;
  beuler->b = NULL;
  if (beuler->s1 != NULL)
    N_VDestroy(beuler->s1);
  beuler->s1 = NULL;
  if (beuler->s2 != NULL)
    N_VDestroy(beuler->s2);
  beuler->s2 = NULL;
  if (beuler->sx != NULL)
    N_VDestroy(beuler->sx);
  beuler->sx = NULL;
  if (beuler->sb != NULL)
    N_VDestroy(beuler->sb);
  beuler->sb = NULL;
}

static void gmres_init(void* context, int N)
{
  linear_beuler_t* beuler = context;
  linear_beuler_reset(beuler);
  beuler->x = N_VNew(beuler->comm, N);
  beuler->b = N_VClone(beuler->x);
  SpgmrMem solver = SpgmrMalloc(beuler->max_kdim, beuler->x);
  beuler->solver = solver;
}

static void gmres_step(void* context, double t1, double t2, double* solution, int N)
{
  linear_beuler_t* beuler = context;
  ASSERT(beuler->solver != NULL);
  SpbcgMem solver = beuler->solver;
  double* x = N_VGetArrayPointer(beuler->x);
  for (int i = 0; i < N; ++i)
    x[i] = solution[i];
  beuler->compute_rhs(beuler->context, t2, beuler->b);
  SpbcgSolve(solver, beuler->context, beuler->x, beuler->b, 
             beuler->precond_type, beuler->delta, beuler->context,
             beuler->sx, beuler->sb, beuler->Ax, beuler->precond,
             &beuler->res_norm, &beuler->nli, &beuler->nps);
}

static void gmres_dtor(void* context)
{
  linear_beuler_t* beuler = context;
  linear_beuler_reset(beuler);
  if ((beuler->context != NULL) && (beuler->dtor != NULL))
    beuler->dtor(beuler->context);
  free(beuler);
}

integrator_t* gmres_linear_beuler_integrator_new(MPI_Comm comm,
                                                 void* context, 
                                                 integrator_Ax_func Ax,
                                                 integrator_compute_rhs_func compute_rhs,
                                                 int precond_type,
                                                 integrator_precond_func precond,
                                                 int max_kdim,
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
  ASSERT(max_kdim > 0);
  ASSERT((gram_schmidt == CLASSICAL_GS) || (gram_schmidt == MODIFIED_GS));
  ASSERT(delta > 0.0);
  ASSERT(max_restarts >= 0);
  linear_beuler_t* beuler = malloc(sizeof(linear_beuler_t));
  beuler->comm = comm;
  beuler->context = context;
  beuler->solver = NULL;
  beuler->Ax = Ax;
  beuler->compute_rhs = compute_rhs;
  beuler->precond_type = precond_type;
  beuler->precond = precond;
  beuler->max_kdim = max_kdim;
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
  linear_beuler_t* beuler = context;
  linear_beuler_reset(beuler);
  beuler->x = N_VNew(beuler->comm, N);
  beuler->b = N_VClone(beuler->x);
  SpbcgMem solver = SpbcgMalloc(beuler->max_kdim, beuler->x);
  beuler->solver = solver;
}

static void bicgstab_step(void* context, double t1, double t2, double* solution, int N)
{
  linear_beuler_t* beuler = context;
  ASSERT(beuler->solver != NULL);
  SpbcgMem solver = beuler->solver;
  beuler->compute_rhs(beuler->context, t2, beuler->b);
  SpbcgSolve(solver, beuler->context, beuler->x, beuler->b, 
             beuler->precond_type, beuler->delta, beuler->context,
             beuler->sx, beuler->sb, beuler->Ax, beuler->precond,
             &beuler->res_norm, &beuler->nli, &beuler->nps);
}

static void bicgstab_dtor(void* context)
{
  linear_beuler_t* beuler = context;
  linear_beuler_reset(beuler);
  if ((beuler->context != NULL) && (beuler->dtor != NULL))
    beuler->dtor(beuler->context);
  free(beuler);
}

integrator_t* bicgstab_linear_beuler_integrator_new(MPI_Comm comm,
                                                    void* context, 
                                                    integrator_Ax_func Ax,
                                                    integrator_compute_rhs_func compute_rhs,
                                                    int precond_type,
                                                    integrator_precond_func precond,
                                                    int max_kdim,
                                                    double delta,
                                                    integrator_dtor dtor)
{
  ASSERT(Ax != NULL);
  ASSERT(compute_rhs != NULL);
  ASSERT((precond_type == PREC_NONE) || (precond_type == PREC_LEFT) ||
         (precond_type == PREC_RIGHT) || (precond_type == PREC_BOTH));
  ASSERT((precond != NULL) || (precond_type == PREC_NONE));
  ASSERT(max_kdim > 0);
  ASSERT(delta > 0.0);
  linear_beuler_t* beuler = malloc(sizeof(linear_beuler_t));
  beuler->comm = comm;
  beuler->context = context;
  beuler->solver = NULL;
  beuler->Ax = Ax;
  beuler->compute_rhs = compute_rhs;
  beuler->precond_type = precond_type;
  beuler->precond = precond;
  beuler->max_kdim = max_kdim;
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

