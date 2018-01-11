// Copyright (c) 2012-2018, Jeffrey N. Johnson
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
#include "cmocka.h"
#include "core/polynomial.h"
#include "core/rng.h"

#if POLYMEC_HAVE_DOUBLE_PRECISION
static const real_t machine_precision = 1e-14;
#else
static const real_t machine_precision = 1e-7;
#endif

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

static void test_ctor(void** state, int p)
{
  // Set up coefficients.
  int dim = polynomial_basis_dim(p);
  real_t coeffs[dim];
  for (int i = 0; i < dim; ++i)
    coeffs[i] = 1.0*i;

  // Construct a polynomial about the origin.
  polynomial_t* poly = polynomial_new(p, coeffs, NULL);
  assert_int_equal(dim, polynomial_num_terms(poly));
  point_t origin = {.x = 0.0, .y = 0.0, .z = 0.0};
  assert_true(point_distance(&origin, polynomial_x0(poly)) < machine_precision);
  int pos = 0, x_pow, y_pow, z_pow, index = 0;
  real_t coeff;
  while (polynomial_next(poly, &pos, &coeff, &x_pow, &y_pow, &z_pow))
  {
    assert_true(reals_equal(coeff, coeffs[index]));
    assert_int_equal(x_powers[p][index], x_pow);
    assert_int_equal(y_powers[p][index], y_pow);
    assert_int_equal(z_powers[p][index], z_pow);
    ++index;
  }

  // Shift it by x0.
  point_t x0 = {.x = 1.0, .y = 2.0, .z = 3.0};
  polynomial_shift(poly, &x0);

  // Scale it by 1.
  polynomial_scale(poly, 1.0);

  // Construct a polynomial about the same point x0 and compare it to the 
  // first one.
  polynomial_t* poly1 = polynomial_new(p, coeffs, &x0);
  polynomial_fprintf(poly, stdout);
  assert_true(polynomial_equals(poly, poly1));
  assert_int_equal(dim, polynomial_num_terms(poly));
  assert_true(point_distance(&x0, polynomial_x0(poly)) < machine_precision);
  pos = 0; index = 0;
  while (polynomial_next(poly, &pos, &coeff, &x_pow, &y_pow, &z_pow))
  {
    assert_true(reals_equal(coeff, coeffs[index]));
    assert_int_equal(x_powers[p][index], x_pow);
    assert_int_equal(y_powers[p][index], y_pow);
    assert_int_equal(z_powers[p][index], z_pow);
    int I = polynomial_index(poly, x_pow, y_pow, z_pow);
    assert_int_equal(index, I);
    ++index;
  }

  // Clone the polynomial and do an equality comparison.
  polynomial_t* poly2 = polynomial_clone(poly1);
  assert_true(polynomial_equals(poly1, poly2));
  poly2 = NULL;

  // Create a polynomial representing the x derivative of this one.
  polynomial_t* dpoly = polynomial_derivative(poly, 1, 0, 0);
  assert_true(reals_equal(polynomial_value(dpoly, &x0), 
                          polynomial_deriv_value(poly, 1, 0, 0, &x0)));

  // Create a polynomial representing the y derivative of this one.
  dpoly = polynomial_derivative(poly, 0, 1, 0);
  assert_true(reals_equal(polynomial_value(dpoly, &x0), 
                          polynomial_deriv_value(poly, 0, 1, 0, &x0)));

  // Create a polynomial representing the z derivative of this one.
  dpoly = polynomial_derivative(poly, 0, 0, 1);
  assert_true(reals_equal(polynomial_value(dpoly, &x0), 
                          polynomial_deriv_value(poly, 0, 0, 1, &x0)));

  // Add this polynomial to itself.
  polynomial_add(poly1, 1.0, poly);

  // Compute the product of the original polynomial with itself.
  polynomial_t* poly3 = polynomial_product(poly, poly);

  poly = NULL;
  poly1 = NULL;
  poly2 = NULL;
  poly3 = NULL;
  dpoly = NULL;
}

static void test_p0_ctor(void** state)
{
  test_ctor(state, 0);
}

static void test_p1_ctor(void** state)
{
  test_ctor(state, 1);
}

static void test_p2_ctor(void** state)
{
  test_ctor(state, 2);
}

static void test_p3_ctor(void** state)
{
  test_ctor(state, 3);
}

static void test_p4_ctor(void** state)
{
  test_ctor(state, 4);
}

static int fact(int n)
{
  if (n <= 0)
    return 1;
  else
    return n * fact(n-1);
}

static void test_basis(void** state, int p)
{
  // We test the polynomial with coefficients 1, 2, 3, 4, ... in the standard
  // basis.
  int dim = polynomial_basis_dim(p);
  real_t coeffs[dim];
  for (int i = 0; i < dim; ++i)
    coeffs[i] = 1.0*i;
  polynomial_t* poly = polynomial_new(p, coeffs, NULL);

  // Evaluate the basis of the polynomial at a random point x.
  rng_t* rng = host_rng_new();
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0,
                 .y1 = 0.0, .y2 = 1.0,
                 .z1 = 0.0, .z2 = 1.0};
  point_t x;
  point_randomize(&x, rng, &bbox);

  // Basis itself.
  {
    real_t basis[dim];
    polynomial_compute_basis(p, 0, 0, 0, &x, basis);

    // Now check everything.
    int pos = 0, x_pow, y_pow, z_pow, index = 0;
    real_t coeff;
    while (polynomial_next(poly, &pos, &coeff, &x_pow, &y_pow, &z_pow))
    {
      real_t term = pow(x.x, x_pow) * pow(x.y, y_pow) * pow(x.z, z_pow);
      assert_true(reals_nearly_equal(term, basis[index], machine_precision));
      ++index;
    }
  }

  // Basis derivatives.
  for (int i = 0; i < p; ++i)
  {
    for (int j = 0; j < p; ++j)
    {
      for (int k = 0; k < p; ++k)
      {
        real_t basis_deriv[dim];
        polynomial_compute_basis(p, i, j, k, &x, basis_deriv);

        int pos = 0, x_pow, y_pow, z_pow, index = 0;
        real_t coeff;
        while (polynomial_next(poly, &pos, &coeff, &x_pow, &y_pow, &z_pow))
        {
          real_t x_term = (x_pow >= i) ? pow(x.x, x_pow - i) * fact(x_pow) / fact(x_pow - i) : 0.0;
          real_t y_term = (y_pow >= j) ? pow(x.y, y_pow - j) * fact(y_pow) / fact(y_pow - j) : 0.0;
          real_t z_term = (z_pow >= k) ? pow(x.z, z_pow - k) * fact(z_pow) / fact(z_pow - k) : 0.0;
          real_t term = (i+j+k > p) ? 0.0 : x_term * y_term * z_term;
          assert_true(reals_nearly_equal(term, basis_deriv[index], machine_precision));
          ++index;
        }
      }
    }
  }

  poly = NULL;
}

static void test_p0_basis(void** state)
{
  test_basis(state, 0);
}

static void test_p1_basis(void** state)
{
  test_basis(state, 1);
}

static void test_p2_basis(void** state)
{
  test_basis(state, 2);
}

static void test_p3_basis(void** state)
{
  test_basis(state, 3);
}

static void test_p4_basis(void** state)
{
  test_basis(state, 4);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);

  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_p0_ctor),
    cmocka_unit_test(test_p1_ctor),
    cmocka_unit_test(test_p2_ctor),
    cmocka_unit_test(test_p3_ctor),
    cmocka_unit_test(test_p4_ctor),
    cmocka_unit_test(test_p0_basis),
    cmocka_unit_test(test_p1_basis),
    cmocka_unit_test(test_p2_basis),
    cmocka_unit_test(test_p3_basis),
    cmocka_unit_test(test_p4_basis)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
