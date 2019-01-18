// Copyright (c) 2012-2019, Jeffrey N. Johnson
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
#include "solvers/ark_ode_solver.h"

extern ode_solver_t* functional_ark_diurnal_solver_new(void);
extern ode_solver_t* bj_jfnk_ark_diurnal_solver_new(newton_pc_side_t side);
extern real_t* diurnal_initial_conditions(ode_solver_t* integ);
extern int test_diurnal_step(void** state, ode_solver_t* integ, int max_step);

static void test_functional_ark_diurnal_ctor(void** state)
{
  ode_solver_t* integ = functional_ark_diurnal_solver_new();
  ode_solver_free(integ);
}

static void test_bj_jfnk_ark_diurnal_ctor(void** state)
{
  ode_solver_t* integ = bj_jfnk_ark_diurnal_solver_new(NEWTON_PC_LEFT);
  ode_solver_free(integ);
  integ = bj_jfnk_ark_diurnal_solver_new(NEWTON_PC_RIGHT);
  ode_solver_free(integ);
  integ = bj_jfnk_ark_diurnal_solver_new(NEWTON_PC_BOTH);
  ode_solver_free(integ);
}

//static void test_functional_ark_diurnal_step(void** state)
//{
//  ode_solver_t* integ = functional_ark_diurnal_solver_new();
//  test_diurnal_step(state, integ, 1.0, 5000);
//}

static void test_bj_jfnk_ark_diurnal_step_left(void** state)
{
  ode_solver_t* integ = bj_jfnk_ark_diurnal_solver_new(NEWTON_PC_LEFT);
#if POLYMEC_HAVE_DOUBLE_PRECISION
  int max_steps = 500;
#else
  int max_steps = 313;
#endif
  test_diurnal_step(state, integ, max_steps);
}

//static void test_bj_jfnk_ark_diurnal_step_right(void** state)
//{
//  ode_solver_t* integ = bj_jfnk_ark_diurnal_solver_new(NEWTON_PC_RIGHT);
//#if POLYMEC_HAVE_DOUBLE_PRECISION
//  int max_steps = 500;
//#else
//  int max_steps = 344;
//#endif
//  test_diurnal_step(state, integ, max_steps);
//}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_functional_ark_diurnal_ctor),
    cmocka_unit_test(test_bj_jfnk_ark_diurnal_ctor),
//    cmocka_unit_test(test_functional_ark_diurnal_step), // too stiff!
    cmocka_unit_test(test_bj_jfnk_ark_diurnal_step_left),
//    cmocka_unit_test(test_bj_jfnk_ark_diurnal_step_right) // not working for some reason
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
