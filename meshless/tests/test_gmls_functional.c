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
#include "make_mlpg_lattice.h"
#include "poisson_gmls_functional.h"

void test_gmls_functional_ctor(void** state, int p)
{
  point_cloud_t* points;
  real_t* subdomain_extents;
  make_mlpg_lattice(10, 10, 10, 1.0*p, &points, &subdomain_extents, NULL);
  real_t delta = 0.5;
  gmls_functional_t* poisson = poisson_gmls_functional_new(p, points, subdomain_extents, delta);
  assert_int_equal(1, gmls_functional_num_components(poisson));

  // Clean up.
  gmls_functional_free(poisson);
  point_cloud_free(points);
  polymec_free(subdomain_extents);
}

void test_gmls_functional_ctor_2(void** state)
{
  test_gmls_functional_ctor(state, 2);
}

void test_gmls_functional_ctor_3(void** state)
{
  test_gmls_functional_ctor(state, 3);
}

void test_gmls_functional_ctor_4(void** state)
{
  test_gmls_functional_ctor(state, 4);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_gmls_functional_ctor_2),
    unit_test(test_gmls_functional_ctor_3),
    unit_test(test_gmls_functional_ctor_4)
  };
  return run_tests(tests);
}
