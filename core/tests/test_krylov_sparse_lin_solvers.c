#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmockery.h"
#include "core/krylov_sparse_lin_solvers.h"

typedef struct 
{
  double dx; // Grid spacing.
  double x1, x2; // Dirichlet BCs for the left and right ends.
} laplacian_test_t;

// This function computes the Laplacian operator applied to the vector x, with 
// Dirichlet boundary conditions at the ends of a 1D domain.
static int compute_Lx(void* context, N_Vector x, N_Vector Lx)
{
  laplacian_test_t* test = context;
  return 0;
}

// This function solves the preconditioner system M * z = r, where M is the diagonal 
// of the Laplacian operator.
static int jacobi_precond(void* context, N_Vector r, N_Vector z, int precond_type)
{
  laplacian_test_t* test = context;
  double* zdata = NV_DATA(z);
  double* rdata = NV_DATA(r);
  int n = NV_LOCLENGTH(z);
  if (precond_type == PREC_BOTH) 
  {
    // We are doing symmetric preconditioning, so Mii is the square root 
    // of the diagonal of the laplacian, and Mij = 0 for j != i.
    for (int i = 0; i < n; ++i)
    {
      double Mii = sqrt(2.0 / test->dx);
      rdata[i] = zdata[i] / Mii;
    }
  }
  else if (precond_type != PREC_NONE)
  {
    // Mii is simply the diagonal of the laplacian.
    for (int i = 0; i < n; ++i)
    {
      double Mii = 2.0 / test->dx;
      rdata[i] = zdata[i] / Mii;
    }
  }
  return 0;
}

void test_laplace_equation_with_solver(void** state, sparse_lin_solver_t* solver)
{
  // FIXME
}

void test_laplace_equation(void** state)
{
  // Grid spacing, Dirichlet BCs.
  laplacian_test_t test = {.dx = 0.01, .x1 = 0.0, .x2 = 1.0};

  // GMRES solver with no preconditioner.
  MPI_Comm comm = MPI_COMM_WORLD;
  int max_kdim = 5;
  double delta = 1e-8;
  int max_restarts = 20;
  sparse_lin_solver_t* solver;
  solver = gmres_sparse_lin_solver_new(comm, &test, compute_Lx, max_kdim, MODIFIED_GS,
                                       PREC_NONE, NULL, delta, max_restarts, NULL);
  test_laplace_equation_with_solver(state, solver);
  sparse_lin_solver_free(solver);

  // GMRES solver with Jacobi (diagonal scaling) preconditioner, applied from the left.
  solver = gmres_sparse_lin_solver_new(comm, &test, compute_Lx, max_kdim, MODIFIED_GS,
                                       PREC_LEFT, jacobi_precond, delta, max_restarts, NULL);
  test_laplace_equation_with_solver(state, solver);
  sparse_lin_solver_free(solver);

  // BiCGStab solver with no preconditioner.
  solver = bicgstab_sparse_lin_solver_new(comm, &test, compute_Lx, max_kdim, 
                                          PREC_NONE, NULL, delta, NULL);
  test_laplace_equation_with_solver(state, solver);
  sparse_lin_solver_free(solver);

  // BiCGStab solver with Jacobi (diagonal scaling) preconditioner, applied from the left.
  solver = bicgstab_sparse_lin_solver_new(comm, &test, compute_Lx, max_kdim, 
                                          PREC_LEFT, jacobi_precond, delta, NULL);
  test_laplace_equation_with_solver(state, solver);
  sparse_lin_solver_free(solver);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_laplace_equation)
  };
  return run_tests(tests);
}
