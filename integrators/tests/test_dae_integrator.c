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

extern dae_integrator_t* bj_pc_jfnk_heat2d_integrator_new(void);
extern dae_integrator_t* ink_heat2d_integrator_new(krylov_factory_t* factory);
extern void heat2d_set_initial_conditions(dae_integrator_t* integ, real_t** u, real_t** u_dot);
extern krylov_factory_t* create_petsc_krylov_factory(void);
extern krylov_factory_t* create_hypre_krylov_factory(void);

static void test_bj_pc_jfnk_heat2d_ctor(void** state)
{
  dae_integrator_t* integ = bj_pc_jfnk_heat2d_integrator_new();
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
  polymec_free(u);
  polymec_free(udot);
}

static void test_bj_pc_jfnk_heat2d_step(void** state)
{
  dae_integrator_t* integ = bj_pc_jfnk_heat2d_integrator_new();
  test_heat2d_step(state, integ);
}

static void test_ink_heat2d_ctor(void** state, krylov_factory_t* factory)
{
  dae_integrator_t* integ = ink_heat2d_integrator_new(factory);
  dae_integrator_free(integ);
}

static void test_petsc_ink_heat2d_ctor(void** state)
{
  krylov_factory_t* petsc = create_petsc_krylov_factory();
  if (petsc != NULL)
    test_ink_heat2d_ctor(state, petsc);
}

static void test_hypre_ink_heat2d_ctor(void** state)
{
  krylov_factory_t* hypre = create_hypre_krylov_factory();
  if (hypre != NULL)
    test_ink_heat2d_ctor(state, hypre);
}

static void test_ink_heat2d_step(void** state, krylov_factory_t* factory)
{
  dae_integrator_t* integ = ink_heat2d_integrator_new(factory);
  test_heat2d_step(state, integ);
}

static void test_petsc_ink_heat2d_step(void** state)
{
  krylov_factory_t* petsc = create_petsc_krylov_factory();
  if (petsc != NULL)
    test_ink_heat2d_step(state, petsc);
}

static void test_hypre_ink_heat2d_step(void** state)
{
  krylov_factory_t* hypre = create_hypre_krylov_factory();
  if (hypre != NULL)
    test_ink_heat2d_step(state, hypre);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_bj_pc_jfnk_heat2d_ctor),
    cmocka_unit_test(test_petsc_ink_heat2d_ctor),
    cmocka_unit_test(test_hypre_ink_heat2d_ctor),
    cmocka_unit_test(test_bj_pc_jfnk_heat2d_step),
    cmocka_unit_test(test_petsc_ink_heat2d_step),
    cmocka_unit_test(test_hypre_ink_heat2d_step)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
