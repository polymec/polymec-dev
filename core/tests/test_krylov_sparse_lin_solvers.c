#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmockery.h"
#include "core/krylov_sparse_lin_solvers.h"

typedef struct 
{
  MPI_Comm comm;
  double dx; // Grid spacing.
  double x1, x2; // Dirichlet boundary values.
} laplacian_test_t;

// This function computes the Laplacian operator applied to the vector x, with 
// Dirichlet boundary conditions at the ends of a 1D domain.
static int compute_Lx(void* context, N_Vector x, N_Vector Lx)
{
  laplacian_test_t* test = context;
  double* xdata = NV_DATA(x);
  double* Lxdata = NV_DATA(Lx);
  int n = NV_LOCLENGTH(x);
  double dx2 = test->dx * test->dx;

  // Apply the general stencil to the interior.
  for (int i = 1; i < n-1; ++i)
    Lxdata[i] = (xdata[i-1] - 2.0*xdata[i] + xdata[i+1]) / dx2;

  // Apply boundary conditions, parallel-related and otherwise.
#if USE_MPI
  int N, rank;
  MPI_Comm_size(test->comm, &N);
  MPI_Comm_rank(test->comm, &rank);
  int tag = 0;
  double xleft, xright;
  if (rank == 0)
  {
    MPI_Comm_send(&xdata[n-1], rank+1, MPI_DOUBLE, 1, tag, test->comm);
    MPI_Comm_recv(&xright, rank+1, MPI_DOUBLE, 1, tag, test->comm);
    Lxdata[0] = xdata[0];
    Lxdata[n-1] = (xdata[i-1] - 2.0*xdata[i] + xright) / dx2;
  }
  else if (rank == N-1)
  {
    MPI_Comm_send(&xdata[0], rank-1, MPI_DOUBLE, 1, tag, test->comm);
    MPI_Comm_recv(&xleft, rank-1, MPI_DOUBLE, 1, tag, test->comm);
    Lxdata[0] = (xleft - 2.0*xdata[i] + xdata[i+1]) / dx2;
    Lxdata[n-1] = xdata[n-1];
  }
  else
  {
    MPI_Comm_recv(&xleft, rank-1, MPI_DOUBLE, 1, tag, test->comm);
    MPI_Comm_recv(&xright, rank+1, MPI_DOUBLE, 1, tag, test->comm);
    Lxdata[0] = (xleft - 2.0*xdata[i] + xdata[i+1]) / test->dx;
    Lxdata[n-1] = (xdata[i-1] - 2.0*xdata[i] + xright) / dx2;
  }
#else
  Lxdata[0] = xdata[0];
  Lxdata[n-1] = xdata[n-1];
#endif

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

void test_laplace_equation_with_solver(void** state, sparse_lin_solver_t* solver, laplacian_test_t* test)
{
  int N = 8;
  N_Vector x = N_VNew(MPI_COMM_WORLD, N);
  N_Vector b = N_VNew(MPI_COMM_WORLD, N);
  double* xdata = NV_DATA(x);
  double* bdata = NV_DATA(b);
  for (int i = 0; i < N; ++i)
    xdata[i] = bdata[i] = 0.0;
#if USE_MPI
  int nproc, rank;
  MPI_Comm_size(test->comm, &nproc);
  MPI_Comm_rank(test->comm, &rank);
  if (rank == 0)
    bdata[0] = test->x1;
  else if (rank == nproc-1)
    bdata[N-1] = test->x2;
#else
  bdata[0] = test->x1;
  bdata[N-1] = test->x2;
#endif
  sparse_lin_solver_outcome_t outcome = sparse_lin_solver_solve(solver, x, b);
//  assert_true(outcome == SPARSE_LIN_SOLVER_CONVERGED);
  double norm;
  int nli, nps;
  sparse_lin_solver_get_info(solver, &norm, &nli, &nps);
  printf("res norm = %g, nli = %d\n", norm, nli);
  printf("x = ");
  for (int i = 0; i < N; ++i)
    printf("%g ", xdata[i]);
  printf("\n");
  printf("b = ");
  for (int i = 0; i < N; ++i)
    printf("%g ", bdata[i]);
  printf("\n");
  N_VDestroy(x);
  N_VDestroy(b);
}

void test_laplace_equation(void** state)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  int nproc;
  MPI_Comm_size(comm, &nproc);
  int N = nproc * 64;

  // Communicator, grid spacing.
  laplacian_test_t test = {.comm = MPI_COMM_WORLD, .dx = 1.0/N, .x1 = 0.0, .x2 = 1.0 };

  // GMRES solver with no preconditioner.
  int max_kdim = 5;
  double delta = 1e-8;
  int max_restarts = 20;
  sparse_lin_solver_t* solver;
  solver = gmres_sparse_lin_solver_new(comm, &test, compute_Lx, max_kdim, MODIFIED_GS,
                                       PREC_NONE, NULL, delta, max_restarts, NULL);
  test_laplace_equation_with_solver(state, solver, &test);
  sparse_lin_solver_free(solver);

  // GMRES solver with Jacobi (diagonal scaling) preconditioner, applied from the left.
  solver = gmres_sparse_lin_solver_new(comm, &test, compute_Lx, max_kdim, MODIFIED_GS,
                                       PREC_LEFT, jacobi_precond, delta, max_restarts, NULL);
  test_laplace_equation_with_solver(state, solver, &test);
  sparse_lin_solver_free(solver);

  // BiCGStab solver with no preconditioner.
  solver = bicgstab_sparse_lin_solver_new(comm, &test, compute_Lx, max_kdim, 
                                          PREC_NONE, NULL, delta, NULL);
  test_laplace_equation_with_solver(state, solver, &test);
  sparse_lin_solver_free(solver);

  // BiCGStab solver with Jacobi (diagonal scaling) preconditioner, applied from the left.
  solver = bicgstab_sparse_lin_solver_new(comm, &test, compute_Lx, max_kdim, 
                                          PREC_LEFT, jacobi_precond, delta, NULL);
  test_laplace_equation_with_solver(state, solver, &test);
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
