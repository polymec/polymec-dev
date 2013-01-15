#include "integrators/linear_beuler_integrator.h"
#include "core/krylov_sparse_lin_solvers.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
  MPI_Comm comm;
  void* context;
  integrator_compute_Ax_func compute_Ax;
  integrator_compute_rhs_func compute_rhs;
  integrator_apply_bcs_func apply_bcs;
  sparse_lin_solver_t* solver;
  integrator_dtor dtor;
  N_Vector x, Ax, b;
  double dt, t2;
} linear_beuler_t;

static void beuler_reset(linear_beuler_t* beuler)
{
  if (beuler->x != NULL)
  {
    N_VDestroy(beuler->x);
    N_VDestroy(beuler->Ax);
    N_VDestroy(beuler->b);
  }
  beuler->x = NULL;
  beuler->Ax = NULL;
  beuler->b = NULL;
}

static void beuler_init(void* context, int N)
{
  linear_beuler_t* beuler = context;
  beuler_reset(beuler);
  beuler->x = N_VNew(beuler->comm, N);
  beuler->Ax = N_VClone(beuler->x);
  beuler->b = N_VClone(beuler->x);
}

static void beuler_step(void* context, double t1, double t2, double* solution, int N)
{
  linear_beuler_t* beuler = context;
  beuler->t2 = t2;
  beuler->dt = t2 - t1;

  // Compute the right hand side at t2.
  beuler->compute_rhs(beuler->context, t2, beuler->b);

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

  // Apply the operator A to the solution vector x.
  if (beuler->compute_Ax(beuler->context, x, beuler->Ax) != 0)
    return -1;

  // Now scale by -dt and add the identity matrix times x.
  N_VLinearSum(1.0, x, -beuler->dt, beuler->Ax, Ax);

  // Finally, apply boundary conditions to Ax.
  if (beuler->apply_bcs != NULL)
    beuler->apply_bcs(beuler->context, beuler->t2, Ax);
  return 0;
}

integrator_t* linear_beuler_integrator_new(MPI_Comm comm,
                                           void* context, 
                                           integrator_compute_Ax_func compute_Ax,
                                           integrator_compute_rhs_func compute_rhs,
                                           integrator_apply_bcs_func apply_bcs,
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
  beuler->apply_bcs = apply_bcs;
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

