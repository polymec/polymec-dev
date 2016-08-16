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
#include "model/physical_constants.h"

#define SI_UNITS 1
#define CGS_UNITS 2

#define BOLTZMANN_CONSTANT 1

void test_define_units_system(void** state)
{
  define_units_system(SI_UNITS);
  define_units_system(CGS_UNITS);
  use_units_system(SI_UNITS);
}

void test_set_physical_constant(void** state)
{
  set_physical_constant(BOLTZMANN_CONSTANT, 1.3807e-23, SI_UNITS);
  set_physical_constant(BOLTZMANN_CONSTANT, 1.3807e-16, CGS_UNITS);
  assert_true(reals_equal(physical_constant_in_units(BOLTZMANN_CONSTANT, CGS_UNITS), 1.3807e-16));

  // Recall that we're still set to the SI units system from the above test.
  assert_true(reals_equal(physical_constant(BOLTZMANN_CONSTANT), 1.3807e-23));
  use_units_system(CGS_UNITS);
  assert_true(reals_equal(physical_constant(BOLTZMANN_CONSTANT), 1.3807e-16));
}

void test_conversion_factor(void** state)
{
  real_t kb1 = physical_constant_in_units(BOLTZMANN_CONSTANT, SI_UNITS);
  real_t kb2 = physical_constant_in_units(BOLTZMANN_CONSTANT, CGS_UNITS);
  real_t factor = physical_constant_conversion_factor(BOLTZMANN_CONSTANT, SI_UNITS, CGS_UNITS);
  assert_true(reals_equal(factor, kb2/kb1));
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_define_units_system),
    cmocka_unit_test(test_set_physical_constant),
    cmocka_unit_test(test_conversion_factor)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
