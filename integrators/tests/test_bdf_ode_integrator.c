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
#include "integrators/bdf_ode_integrator.h"


extern ode_integrator_t* block_jacobi_precond_bdf_diurnal_integrator_new();
extern ode_integrator_t* lu_precond_bdf_diurnal_integrator_new();
extern real_t* diurnal_initial_conditions(ode_integrator_t* integ);

static void test_block_jacobi_precond_diurnal_ctor(void** state)
{
  ode_integrator_t* integ = block_jacobi_precond_bdf_diurnal_integrator_new(NEWTON_PC_LEFT);
  ode_integrator_free(integ);
  integ = block_jacobi_precond_bdf_diurnal_integrator_new(NEWTON_PC_RIGHT);
  ode_integrator_free(integ);
  integ = block_jacobi_precond_bdf_diurnal_integrator_new(NEWTON_PC_BOTH);
  ode_integrator_free(integ);
}

static int test_diurnal_step(void** state, ode_integrator_t* integ, int max_steps)
{
  // Set up the problem.
#if POLYMEC_HAVE_DOUBLE_PRECISION
  bdf_ode_integrator_set_tolerances(integ, 1e-5, 1e-3);
#else
  bdf_ode_integrator_set_tolerances(integ, 1e-4, 1e-2);
#endif
  real_t* u = diurnal_initial_conditions(integ);

  // Integrate it out to t = 86400 s (24 hours).
  real_t t = 0.0;
  int step = 0;
  while (t < 86400.0)
  {
    log_detail("Step %d: t = %g", step, t);
    bool integrated = ode_integrator_step(integ, 7200.0, &t, u);
    assert_true(integrated);
    ++step;

    if (step >= max_steps)
      break;
  }
//printf("u = [");
//for (int i = 0; i < 200; ++i)
//printf("%g ", u[i]);
//printf("]\n");
  printf("Final time: %g\n", t);
  bdf_ode_integrator_diagnostics_t diags;
  bdf_ode_integrator_get_diagnostics(integ, &diags);
  bdf_ode_integrator_diagnostics_fprintf(&diags, stdout);
  assert_true(step < max_steps);

  ode_integrator_free(integ);
  free(u);
  return step;
}

static void test_block_jacobi_precond_diurnal_step_left(void** state)
{
  ode_integrator_t* integ = block_jacobi_precond_bdf_diurnal_integrator_new(NEWTON_PC_LEFT);
#if POLYMEC_HAVE_DOUBLE_PRECISION
  int max_steps = 650;
#else
  int max_steps = 426;
#endif
  test_diurnal_step(state, integ, max_steps);
}

static void test_block_jacobi_precond_diurnal_step_right(void** state)
{
  ode_integrator_t* integ = block_jacobi_precond_bdf_diurnal_integrator_new(NEWTON_PC_RIGHT);
#if POLYMEC_HAVE_DOUBLE_PRECISION
  int max_steps = 500;
#else
  int max_steps = 381;
#endif
  test_diurnal_step(state, integ, max_steps);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_block_jacobi_precond_diurnal_ctor),
    cmocka_unit_test(test_block_jacobi_precond_diurnal_step_left),
    cmocka_unit_test(test_block_jacobi_precond_diurnal_step_right),
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
