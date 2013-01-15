#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmockery.h"
#include "core/krylov_sparse_lin_solvers.h"

void test_laplace_equation_with_solver(void** state, sparse_lin_solver_t* solver)
{
  // FIXME
}

void test_laplace_equation(void** state)
{
  sparse_lin_solver_t* gmres = NULL; //gmres_sparse_lin_solver_new();
  test_laplace_equation_with_solver(state, gmres);
  sparse_lin_solver_free(gmres);
  sparse_lin_solver_t* bicgstab = NULL; //bicgstab_sparse_lin_solver_new();
  test_laplace_equation_with_solver(state, bicgstab);
  sparse_lin_solver_free(bicgstab);
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
