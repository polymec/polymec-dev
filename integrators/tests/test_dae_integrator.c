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
#include "integrators/dae_integrator.h"

extern dae_integrator_t* block_jacobi_precond_heat2d_integrator_new();
extern dae_integrator_t* lu_precond_heat2d_integrator_new();
extern void heat2d_set_initial_conditions(dae_integrator_t* integ, real_t** u, real_t** u_dot);

static void test_block_jacobi_precond_heat2d_ctor(void** state)
{
  dae_integrator_t* integ = block_jacobi_precond_heat2d_integrator_new();
  dae_integrator_free(integ);
}

static void test_lu_precond_heat2d_ctor(void** state)
{
  dae_integrator_t* integ = lu_precond_heat2d_integrator_new();
  dae_integrator_free(integ);
}

static void test_heat2d_step(void** state, dae_integrator_t* integ)
{
  // Set up the problem.
  dae_integrator_set_tolerances(integ, 1e-5, 1e-3);
  real_t *u, *udot;
  heat2d_set_initial_conditions(integ, &u, &udot);

  // Integrate it.
  real_t t = 0.0, tf = 10.24;
  dae_integrator_set_stop_time(integ, tf);
  while (t < tf)
  {
    bool integrated = dae_integrator_step(integ, tf, &t, u, udot);
//    preconditioner_matrix_fprintf(dae_integrator_preconditioner_matrix(integ), stdout);
    assert_true(integrated);
  }
printf("u = [");
for (int i = 0; i < 100; ++i)
printf("%g ", u[i]);
printf("]\n");
  dae_integrator_diagnostics_t diags;
  dae_integrator_get_diagnostics(integ, &diags);
  dae_integrator_diagnostics_fprintf(&diags, stdout);

  dae_integrator_free(integ);
  free(u);
  free(udot);
}

static void test_block_jacobi_precond_heat2d_step(void** state)
{
  dae_integrator_t* integ = block_jacobi_precond_heat2d_integrator_new();
  test_heat2d_step(state, integ);
}

#if POLYMEC_HAVE_DOUBLE_PRECISION
static void test_lu_precond_heat2d_step(void** state)
{
  dae_integrator_t* integ = lu_precond_heat2d_integrator_new();
  test_heat2d_step(state, integ);
}
#endif

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_block_jacobi_precond_heat2d_ctor),
    cmocka_unit_test(test_lu_precond_heat2d_ctor),
    cmocka_unit_test(test_block_jacobi_precond_heat2d_step),
#if POLYMEC_HAVE_DOUBLE_PRECISION
    cmocka_unit_test(test_lu_precond_heat2d_step)
#endif
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
