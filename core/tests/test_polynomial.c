// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <stdlib.h>
#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmockery.h"
#include "core/polynomial.h"

// Powers of x, y, and z for the various polynomial degrees.
static const int x_powers[5][35] = {{0},
                                    {0, 1, 0, 0},
                                    {0, 1, 0, 0, 2, 1, 1, 0, 0, 0},
                                    {0, 1, 0, 0, 2, 1, 1, 0, 0, 0, 3, 2, 2, 1, 1, 1, 0, 0, 0, 0},
                                    {0, 1, 0, 0, 2, 1, 1, 0, 0, 0, 3, 2, 2, 1, 1, 1, 0, 0, 0, 0, 4, 3, 3, 2, 2, 2, 1, 1, 1, 1, 0, 0, 0, 0, 0}};
static const int y_powers[5][35] = {{0},
                                    {0, 0, 1, 0},
                                    {0, 0, 1, 0, 0, 1, 0, 2, 1, 0},
                                    {0, 0, 1, 0, 0, 1, 0, 2, 1, 0, 0, 1, 0, 2, 1, 0, 3, 2, 1, 0},
                                    {0, 0, 1, 0, 0, 1, 0, 2, 1, 0, 0, 1, 0, 2, 1, 0, 3, 2, 1, 0, 0, 1, 0, 2, 1, 0, 3, 2, 1, 0, 4, 3, 2, 1, 0}};
static const int z_powers[5][35] = {{0},
                                    {0, 0, 0, 1},
                                    {0, 0, 0, 1, 0, 0, 1, 0, 1, 2},
                                    {0, 0, 0, 1, 0, 0, 1, 0, 1, 2, 0, 0, 1, 0, 1, 2, 0, 1, 2, 3},
                                    {0, 0, 0, 1, 0, 0, 1, 0, 1, 2, 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4}};

void test_ctor(void** state, int p)
{
  // Set up coefficients.
  int dim = polynomial_basis_dim(p);
  double coeffs[dim];
  for (int i = 0; i < dim; ++i)
    coeffs[i] = 1.0*i;

  // Construct a polynomial about the origin.
  polynomial_t* poly = polynomial_new(p, coeffs, NULL);
  assert_int_equal(dim, polynomial_num_terms(poly));
  point_t origin = {.x = 0.0, .y = 0.0, .z = 0.0};
  assert_true(point_distance(&origin, polynomial_x0(poly)) < 1e-14);
  int pos = 0, x_pow, y_pow, z_pow, index = 0;
  double coeff;
  while (polynomial_next(poly, &pos, &coeff, &x_pow, &y_pow, &z_pow))
  {
    assert_true(coeff == coeffs[index]);
    assert_int_equal(x_powers[p][index], x_pow);
    assert_int_equal(y_powers[p][index], y_pow);
    assert_int_equal(z_powers[p][index], z_pow);
    ++index;
  }

  // Construct a polynomial about a point x0.
  point_t x0 = {.x = 1.0, .y = 2.0, .z = 3.0};
  poly = polynomial_new(p, coeffs, &x0);
  assert_int_equal(dim, polynomial_num_terms(poly));
  assert_true(point_distance(&x0, polynomial_x0(poly)) < 1e-14);
  pos = 0, index = 0;
  while (polynomial_next(poly, &pos, &coeff, &x_pow, &y_pow, &z_pow))
  {
    assert_true(coeff == coeffs[index]);
    assert_int_equal(x_powers[p][index], x_pow);
    assert_int_equal(y_powers[p][index], y_pow);
    assert_int_equal(z_powers[p][index], z_pow);
    ++index;
  }

  poly = NULL;
}

void test_p0_ctor(void** state)
{
  test_ctor(state, 0);
}

void test_p1_ctor(void** state)
{
  test_ctor(state, 1);
}

void test_p2_ctor(void** state)
{
  test_ctor(state, 2);
}

void test_p3_ctor(void** state)
{
  test_ctor(state, 3);
}

void test_p4_ctor(void** state)
{
  test_ctor(state, 3);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);

  const UnitTest tests[] = 
  {
    unit_test(test_p0_ctor),
    unit_test(test_p1_ctor),
    unit_test(test_p2_ctor),
    unit_test(test_p3_ctor),
    unit_test(test_p4_ctor)
  };
  return run_tests(tests);
}
