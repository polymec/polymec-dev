// Copyright (c) 2012-2017, Jeffrey N. Johnson
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
#include "integrators/ark_ode_integrator.h"

extern ode_integrator_t* ink_ark_diurnal_integrator_new(krylov_factory_t* factory);
extern real_t* diurnal_initial_conditions(ode_integrator_t* integ);
extern krylov_factory_t* create_petsc_krylov_factory(void);
extern krylov_factory_t* create_hypre_krylov_factory(void);
extern int test_diurnal_step(void** state, ode_integrator_t* integ, int max_steps);

static void test_ink_ark_diurnal_ctor(void** state, krylov_factory_t* factory)
{
  ode_integrator_t* integ = ink_ark_diurnal_integrator_new(factory);
  ode_integrator_free(integ);
}

static void test_petsc_ink_ark_diurnal_ctor(void** state)
{
  krylov_factory_t* petsc = create_petsc_krylov_factory();
  if (petsc != NULL)
    test_ink_ark_diurnal_ctor(state, petsc);
}

static void test_hypre_ink_ark_diurnal_ctor(void** state)
{
  krylov_factory_t* hypre = create_hypre_krylov_factory();
  if (hypre != NULL)
    test_ink_ark_diurnal_ctor(state, hypre);
}

static void test_ink_ark_diurnal_step(void** state, krylov_factory_t* factory)
{
  ode_integrator_t* integ = ink_ark_diurnal_integrator_new(factory);
#if POLYMEC_HAVE_DOUBLE_PRECISION
  int max_steps = 550;
#else
  int max_steps = 450;
#endif
  test_diurnal_step(state, integ, max_steps);
}

static void test_petsc_ink_ark_diurnal_step(void** state)
{
  krylov_factory_t* petsc = create_petsc_krylov_factory();
  if (petsc != NULL)
    test_ink_ark_diurnal_step(state, petsc);
}

static void test_hypre_ink_ark_diurnal_step(void** state)
{
  krylov_factory_t* hypre = create_hypre_krylov_factory();
  if (hypre != NULL)
    test_ink_ark_diurnal_step(state, hypre);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_petsc_ink_ark_diurnal_ctor),
    cmocka_unit_test(test_hypre_ink_ark_diurnal_ctor),
    cmocka_unit_test(test_petsc_ink_ark_diurnal_step),
    cmocka_unit_test(test_hypre_ink_ark_diurnal_step)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
