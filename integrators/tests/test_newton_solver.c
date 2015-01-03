// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>

#include "cmockery.h"
#include "core/polymec.h"
#include "core/norms.h"
#include "integrators/newton_solver.h"

extern newton_solver_t* block_jacobi_precond_foodweb_solver_new();
extern newton_solver_t* lu_precond_foodweb_solver_new();
extern newton_solver_t* ilu_precond_foodweb_solver_new();
extern real_t* foodweb_initial_conditions();

void test_block_jacobi_precond_foodweb_ctor(void** state)
{
  newton_solver_t* integ = block_jacobi_precond_foodweb_solver_new();
  assert_true(strcmp(newton_solver_name(integ), "Food web") == 0);
  newton_solver_free(integ);
}

void test_lu_precond_foodweb_ctor(void** state)
{
  newton_solver_t* integ = lu_precond_foodweb_solver_new();
  assert_true(strcmp(newton_solver_name(integ), "Food web") == 0);
  newton_solver_free(integ);
}

void test_ilu_precond_foodweb_ctor(void** state)
{
  newton_solver_t* integ = ilu_precond_foodweb_solver_new();
  assert_true(strcmp(newton_solver_name(integ), "Food web") == 0);
  newton_solver_free(integ);
}

void test_foodweb_solve(void** state, newton_solver_t* integ)
{
  // Set up the problem.
  newton_solver_set_tolerances(integ, 1e-7, 1e-13);
  real_t* cc = foodweb_initial_conditions();

  // Solve it.
  int num_iters;
  bool solved = newton_solver_solve(integ, 0.0, cc, &num_iters);
  if (!solved)
  {
    newton_solver_diagnostics_t diagnostics;
    newton_solver_get_diagnostics(integ, &diagnostics);
    newton_solver_diagnostics_fprintf(&diagnostics, stdout);
    preconditioner_fprintf(newton_solver_preconditioner(integ), stdout);
  }
  assert_true(solved);
  log_info("num iterations = %d\n", num_iters);
  assert_true(num_iters < 10);

  // Evaluate the 2-norm of the residual.
  int num_eq = 6*8*8;
  real_t F[num_eq];
  newton_solver_eval_residual(integ, 0.0, cc, F);
  real_t L2 = l2_norm(F, num_eq);
  log_info("||F||_L2 = %g\n", L2);
  assert_true(L2 < sqrt(1e-7));

  newton_solver_free(integ);
  free(cc);
}

void test_block_jacobi_precond_foodweb_solve(void** state)
{
  // Set up the problem.
  newton_solver_t* integ = block_jacobi_precond_foodweb_solver_new();
  test_foodweb_solve(state, integ);
}

void test_lu_precond_foodweb_solve(void** state)
{
  // Set up the problem.
  newton_solver_t* integ = lu_precond_foodweb_solver_new();
  test_foodweb_solve(state, integ);
}

void test_ilu_precond_foodweb_solve(void** state)
{
  // Set up the problem.
  newton_solver_t* integ = ilu_precond_foodweb_solver_new();
  test_foodweb_solve(state, integ);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_block_jacobi_precond_foodweb_ctor),
    unit_test(test_lu_precond_foodweb_ctor),
    unit_test(test_ilu_precond_foodweb_ctor),
    unit_test(test_block_jacobi_precond_foodweb_solve),
    unit_test(test_lu_precond_foodweb_solve),
//    unit_test(test_ilu_precond_foodweb_solve),
  };
  return run_tests(tests);
}
