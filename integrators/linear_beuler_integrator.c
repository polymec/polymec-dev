#include "integrators/linear_beuler_integrator.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
  MPI_Comm comm;
  void* context;
  integrator_compute_Ax_func compute_Ax;
  integrator_compute_rhs_func compute_rhs;
  sparse_lin_solver_t* solver;
  integrator_dtor dtor;
  N_Vector x, b;
} linear_beuler_t;

static void beuler_reset(linear_beuler_t* beuler)
{
  if (beuler->x != NULL)
  {
    N_VDestroy(beuler->x);
    N_VDestroy(beuler->b);
  }
  beuler->x = NULL;
  beuler->b = NULL;
}

static void beuler_init(void* context, int N)
{
  linear_beuler_t* beuler = context;
  beuler_reset(beuler);
  beuler->x = N_VNew(beuler->comm, N);
  beuler->b = N_VClone(beuler->x);
}

static void beuler_step(void* context, double t1, double t2, double* solution, int N)
{
  linear_beuler_t* beuler = context;
  // FIXME: This probably doesn't work.
  double* x = N_VGetArrayPointer(beuler->x);
  for (int i = 0; i < N; ++i)
    x[i] = solution[i];
  sparse_lin_solver_solve(beuler->solver, beuler->x, beuler->b);
  // FIXME
}

static void beuler_dtor(void* context)
{
  linear_beuler_t* beuler = context;
  beuler_reset(beuler);
  if ((beuler->context != NULL) && (beuler->dtor != NULL))
    beuler->dtor(beuler->context);
  sparse_lin_solver_free(beuler->solver);
  free(beuler);
}

static int beuler_compute_Ax(void* context, N_Vector x, N_Vector Ax)
{
  linear_beuler_t* beuler = context;
  // FIXME
  return 0;
}

integrator_t* linear_beuler_integrator_new(MPI_Comm comm,
                                           void* context, 
                                           integrator_compute_Ax_func compute_Ax,
                                           integrator_compute_rhs_func compute_rhs,
                                           integrator_dtor dtor)
{
  ASSERT(compute_Ax != NULL);
  ASSERT(compute_rhs != NULL);
  ASSERT(solver != NULL);
  linear_beuler_t* beuler = malloc(sizeof(linear_beuler_t));
  beuler->comm = comm;
  beuler->context = context;
  beuler->compute_Ax = compute_Ax;
  beuler->compute_rhs = compute_rhs;
  beuler->solver = NULL; // FIXME
  beuler->dtor = dtor;
  integrator_vtable vtable = { .init = beuler_init, .step = beuler_step, .dtor = beuler_dtor };
  return integrator_new("Backward Euler", beuler, vtable, 1, INTEGRATOR_IMPLICIT);
}

sparse_lin_solver_t* linear_beuler_integrator_get_solver(integrator_t* integrator)
{
  linear_beuler_t* beuler = integrator_context(integrator);
  return beuler->solver;
}

#ifdef __cplusplus
}
#endif

