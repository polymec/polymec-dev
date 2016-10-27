// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>

#include "cmocka.h"
#include "core/polymec.h"
#include "core/norms.h"
#include "integrators/newton_solver.h"

extern newton_solver_t* block_jacobi_precond_foodweb_solver_new(void);
extern newton_solver_t* ink_foodweb_solver_new(krylov_factory_t* factory);
extern krylov_factory_t* create_petsc_krylov_factory(void);
extern krylov_factory_t* create_hypre_krylov_factory(void);
extern real_t* foodweb_initial_conditions();

static void test_block_jacobi_precond_foodweb_ctor(void** state)
{
  newton_solver_t* newton = block_jacobi_precond_foodweb_solver_new();
  newton_solver_free(newton);
}

static void test_lis_ink_foodweb_ctor(void** state)
{
  krylov_factory_t* factory = lis_krylov_factory();
  newton_solver_t* newton = ink_foodweb_solver_new(factory);
  newton_solver_free(newton);
}

static void test_petsc_ink_foodweb_ctor(void** state)
{
  krylov_factory_t* factory = create_petsc_krylov_factory();
  if (factory != NULL)
  {
    newton_solver_t* newton = ink_foodweb_solver_new(factory);
    newton_solver_free(newton);
  }
}

static void test_hypre_ink_foodweb_ctor(void** state)
{
  krylov_factory_t* factory = create_hypre_krylov_factory();
  if (factory != NULL)
  {
    newton_solver_t* newton = ink_foodweb_solver_new(factory);
    newton_solver_free(newton);
  }
}

// FIXME: It looks like the food web problem was not tested at single 
// FIXME: precision by the Sundials folks, so we can't use this test 
// FIXME: without double precision.
#if POLYMEC_HAVE_DOUBLE_PRECISION
static void test_foodweb_solve(void** state, newton_solver_t* newton)
{
  // Set up the problem.
  newton_solver_set_tolerances(newton, 1e-7, 1e-13);
  real_t* cc = foodweb_initial_conditions();

  // Solve it.
  int num_iters;
  bool solved = newton_solver_solve(newton, 0.0, cc, &num_iters);
  if (!solved)
  {
    newton_solver_diagnostics_t diagnostics;
    newton_solver_get_diagnostics(newton, &diagnostics);
    newton_solver_diagnostics_fprintf(&diagnostics, stdout);
  }
  assert_true(solved);
  log_info("num iterations = %d\n", num_iters);

  // Evaluate the 2-norm of the residual.
  // Note that the tolerance specifies a bound on the *scaled* norm of F.
  int num_eq = 6*8*8;
  real_t F[num_eq];
  newton_solver_eval_residual(newton, 0.0, cc, F);
  real_t L2 = l2_norm(F, num_eq);
  log_info("||F||_L2 = %g\n", L2);
  assert_true(L2 < 1e-2);

  newton_solver_free(newton);
  polymec_free(cc);
}

static void test_block_jacobi_precond_foodweb_solve(void** state)
{
  // Set up the problem.
  newton_solver_t* newton = block_jacobi_precond_foodweb_solver_new();
  test_foodweb_solve(state, newton);
}

static void test_lis_ink_foodweb_solve(void** state)
{
  krylov_factory_t* factory = lis_krylov_factory();
  newton_solver_t* newton = ink_foodweb_solver_new(factory);
  test_foodweb_solve(state, newton);
}

static void test_petsc_ink_foodweb_solve(void** state)
{
  krylov_factory_t* factory = create_petsc_krylov_factory();
  if (factory != NULL)
  {
    newton_solver_t* newton = ink_foodweb_solver_new(factory);
    test_foodweb_solve(state, newton);
  }
}

static void test_hypre_ink_foodweb_solve(void** state)
{
  krylov_factory_t* factory = create_hypre_krylov_factory();
  if (factory != NULL)
  {
    newton_solver_t* newton = ink_foodweb_solver_new(factory);
    test_foodweb_solve(state, newton);
  }
}

#endif

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_block_jacobi_precond_foodweb_ctor),
    cmocka_unit_test(test_lis_ink_foodweb_ctor),
    cmocka_unit_test(test_petsc_ink_foodweb_ctor),
    cmocka_unit_test(test_hypre_ink_foodweb_ctor),
#if POLYMEC_HAVE_DOUBLE_PRECISION
    cmocka_unit_test(test_block_jacobi_precond_foodweb_solve),
    cmocka_unit_test(test_lis_ink_foodweb_solve),
    cmocka_unit_test(test_petsc_ink_foodweb_solve),
    cmocka_unit_test(test_hypre_ink_foodweb_solve)
#endif
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
