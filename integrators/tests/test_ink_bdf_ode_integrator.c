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

extern ode_integrator_t* ink_bdf_diurnal_integrator_new();
extern real_t* diurnal_initial_conditions(ode_integrator_t* integ);
krylov_factory_t* create_petsc_krylov_factory(void);
krylov_factory_t* create_hypre_krylov_factory(void);

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
    real_t old_t = t;
    bool integrated = ode_integrator_step(integ, 7200.0, &t, u);
    log_detail("Step %d: t = %g, dt = %g", step, t, t - old_t);
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

static void test_ink_bdf_diurnal_ctor(void** state, krylov_factory_t* factory)
{
  ode_integrator_t* integ = ink_bdf_diurnal_integrator_new(factory);
  ode_integrator_free(integ);
}

static void test_lis_ink_bdf_diurnal_ctor(void** state)
{
  krylov_factory_t* lis = lis_krylov_factory();
  test_ink_bdf_diurnal_ctor(state, lis);
}

static void test_petsc_ink_bdf_diurnal_ctor(void** state)
{
  krylov_factory_t* petsc = create_petsc_krylov_factory();
  if (petsc != NULL)
    test_ink_bdf_diurnal_ctor(state, petsc);
}

static void test_hypre_ink_bdf_diurnal_ctor(void** state)
{
  krylov_factory_t* hypre = create_hypre_krylov_factory();
  if (hypre != NULL)
    test_ink_bdf_diurnal_ctor(state, hypre);
}

static void test_ink_bdf_diurnal_step(void** state, krylov_factory_t* factory)
{
  ode_integrator_t* integ = ink_bdf_diurnal_integrator_new(factory);
#if POLYMEC_HAVE_DOUBLE_PRECISION
  int max_steps = 500;
#else
  int max_steps = 381;
#endif
  test_diurnal_step(state, integ, max_steps);
}

static void test_lis_ink_bdf_diurnal_step(void** state)
{
  krylov_factory_t* lis = lis_krylov_factory();
  test_ink_bdf_diurnal_step(state, lis);
}

static void test_petsc_ink_bdf_diurnal_step(void** state)
{
  krylov_factory_t* petsc = create_petsc_krylov_factory();
  if (petsc != NULL)
    test_ink_bdf_diurnal_step(state, petsc);
}

static void test_hypre_ink_bdf_diurnal_step(void** state)
{
  krylov_factory_t* hypre = create_hypre_krylov_factory();
  if (hypre != NULL)
    test_ink_bdf_diurnal_step(state, hypre);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_lis_ink_bdf_diurnal_ctor),
    cmocka_unit_test(test_petsc_ink_bdf_diurnal_ctor),
    cmocka_unit_test(test_hypre_ink_bdf_diurnal_ctor),
    cmocka_unit_test(test_lis_ink_bdf_diurnal_step),
    cmocka_unit_test(test_petsc_ink_bdf_diurnal_step),
    cmocka_unit_test(test_hypre_ink_bdf_diurnal_step)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
