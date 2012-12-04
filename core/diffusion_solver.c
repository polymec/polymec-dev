#include "core/diffusion_solver.h"

#ifdef __cplusplus
extern "C" {
#endif

struct diffusion_solver_t
{
  MPI_Comm comm;

  KSP ksp;
  Mat A;
  Vec b;

  // Flag is set to true if the linear system above is initialized.
  bool initialized;

  // A function for applying boundary conditions.
  diffusion_solver_apply_bcs_func apply_bcs;
  void* context;
  void (*context_dtor)(void*);
};

diffusion_solver_t* diffusion_solver_new(MPI_Comm comm,
                                         diffusion_solver_apply_bcs_func apply_bcs, 
                                         void* context,
                                         void (*context_dtor)(void*))
{
  ASSERT(apply_bcs != NULL);
  diffusion_solver_t* solver = malloc(sizeof(diffusion_solver_t));
  solver->initialized = false;
  solver->comm = comm;
  solver->apply_bcs = apply_bcs;
  solver->context = context;
  solver->context_dtor = context_dtor;
  return solver;
}

void diffusion_solver_free(diffusion_solver_t* solver)
{
  if (solver->initialized)
  {
    KSPDestroy(&solver->ksp);
    MatDestroy(&solver->A);
    VecDestroy(&solver->b);
  }

  if ((solver->context != NULL) && (solver->context_dtor != NULL))
    solver->context_dtor(solver->context);

  free(solver);
}

static void initialize(diffusion_solver_t* solver, 
                       Mat matrix,
                       Vec vector)
{
  if (!solver->initialized)
  {
    KSPCreate(solver->comm, &solver->ksp);
    MatDuplicate(matrix, MAT_COPY_VALUES, &solver->A);
    VecDuplicate(vector, &solver->b);
    solver->initialized = true;
  }
}

static void apply_bcs(diffusion_solver_t* solver, double t)
{
  solver->apply_bcs(solver->context, solver->A, solver->b, t);
}

static void solve(diffusion_solver_t* solver, Vec x)
{
  KSPSetOperators(solver->ksp, solver->A, solver->A, SAME_NONZERO_PATTERN);
  KSPSolve(solver->ksp, solver->b, x);
}

void diffusion_solver_euler(diffusion_solver_t* solver,
                            Mat diffusion_op, 
                            Vec source, 
                            double t1,
                            double t2,
                            Vec sol1, 
                            Vec sol2)
{
  ASSERT(t2 > t1);
  double dt = t2 - t1;

  // Make sure the solver is initialized.
  initialize(solver, diffusion_op, sol1);

  // A = -dt * diffusion_op.
  MatCopy(diffusion_op, solver->A, SAME_NONZERO_PATTERN);
  MatScale(solver->A, -dt);

  // A -> A + I.
  MatShift(solver->A, 1.0);

  // b = sol1 + dt * source.
  VecWAXPY(solver->b, dt, source, sol1);

  // Apply boundary conditions.
  apply_bcs(solver, t2);

  // Solve the linear system.
  solve(solver, sol2);
}

void diffusion_solver_tga(diffusion_solver_t* solver,
                          Mat diffusion_op, 
                          Vec source1, 
                          Vec source2, 
                          double t1,
                          double t2,
                          Vec sol1, 
                          Vec sol2)
{
  ASSERT(t2 > t1);
  double dt = t2 - t1;

  // Make sure the solver is initialized.
  initialize(solver, diffusion_op, sol1);

  // Parameters for the TGA algorithm.
  double a = 2.0 - sqrt(2.0);
  double b = a - 0.5;
  double r1 = (2*a - 1.0) / (a + sqrt(a*a - 4*a + 2.0));
  double r2 = (2*a - 1.0) / (a - sqrt(a*a - 4*a + 2.0));

  // We will do all our work in the new_sol and work vectors.

  // e = [I + (1-a) * dt * D] * old_sol.

  // e -> e + 0.5 * dt * [old_source

  // Do the first solve.
  solve(solver, sol2);

  // Do the second solve.
  solve(solver, sol2);
}

#ifdef __cplusplus
}
#endif

