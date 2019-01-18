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
#include "core/rng.h"
#include "core/lookup1.h"

#define assert_approx_equal(x, y, tol) assert_true(reals_nearly_equal((x), (y), tol))

static void test_zero(void** state, lookup1_interp_t interpolation)
{
  rng_t* rng = host_rng_new();

  int N = 101;
  real_t z[N];
  memset(z, 0, sizeof(real_t) * N);
  lookup1_t* table = lookup1_new(0.0, 1.0, N, z, interpolation);
  assert_true(table != NULL);
  assert_true(lookup1_interpolation(table) == interpolation);
  assert_int_equal(lookup1_num_values(table), N);
  real_t x_min, x_max;
  lookup1_get_bounds(table, &x_min, &x_max);
  assert_true(reals_equal(x_min, 0.0));
  assert_true(reals_equal(x_max, 1.0));
  real_t y_min, y_max;
  lookup1_get_extreme_values(table, &y_min, &y_max);
  assert_true(reals_equal(y_min, 0.0));
  assert_true(reals_equal(y_max, 0.0));
  assert_true(reals_equal(lookup1_value(table, 0.0), 0.0));
  assert_true(reals_equal(lookup1_value(table, 0.5), 0.0));
  assert_true(reals_equal(lookup1_value(table, 1.0), 0.0));
  assert_true(reals_equal(lookup1_value(table, rng_uniform_positive(rng)), 0.0));
  lookup1_free(table);
}

static void test_zero_linear(void** state)
{
  test_zero(state, LOOKUP1_LINEAR);
}

static void test_zero_quadratic(void** state)
{
  test_zero(state, LOOKUP1_QUADRATIC);
}

static void test_line(void** state, lookup1_interp_t interpolation)
{
  rng_t* rng = host_rng_new();

  int N = 101;
  real_t z[N];
  real_t dx = 1.0 / (N-1);
  for (int i = 0; i < N; ++i)
    z[i] = i * dx;
  lookup1_t* table = lookup1_new(0.0, 1.0, N, z, interpolation);
  assert_true(table != NULL);
  assert_int_equal(lookup1_num_values(table), N);
  real_t x_min, x_max;
  lookup1_get_bounds(table, &x_min, &x_max);
  assert_true(reals_equal(x_min, 0.0));
  assert_true(reals_equal(x_max, 1.0));
  real_t y_min, y_max;
  lookup1_get_extreme_values(table, &y_min, &y_max);
  assert_true(reals_equal(y_min, 0.0));
  assert_true(reals_equal(y_max, 1.0));
  assert_true(reals_equal(lookup1_value(table, 0.0), 0.0));
  assert_true(reals_nearly_equal(lookup1_value(table, 0.5), 0.5, 2.0*REAL_EPSILON));
  assert_true(reals_equal(lookup1_value(table, 1.0), 1.0));
  real_t x = rng_uniform_positive(rng);
  assert_approx_equal(lookup1_value(table, x), x, 2.0*REAL_EPSILON);
  lookup1_free(table);
}

static void test_line_linear(void** state)
{
  test_line(state, LOOKUP1_LINEAR);
}

static void test_line_quadratic(void** state)
{
  test_line(state, LOOKUP1_QUADRATIC);
}

static void test_v(void** state, lookup1_interp_t interpolation)
{
  rng_t* rng = host_rng_new();

  int N = 101;
  real_t z[N];
  real_t dx = 1.0 / (N-1);
  for (int i = 0; i < N/2; ++i)
    z[i] = i * dx;
  for (int i = N/2; i < N; ++i)
    z[i] = 0.5 - (i-N/2) * dx;
  lookup1_t* table = lookup1_new(0.0, 1.0, N, z, interpolation);
  assert_true(table != NULL);
  assert_int_equal(lookup1_num_values(table), N);
  assert_true(reals_equal(lookup1_value(table, 0.0), 0.0));
  assert_true(reals_nearly_equal(lookup1_value(table, 0.5), 0.5, 2.0*REAL_EPSILON));
  assert_true(reals_equal(lookup1_value(table, 1.0), 0.0));

  real_t x = rng_uniform_positive(rng);
  // If we are using quadratic interpolation, we have to stay clear of the kink.
  if (interpolation == LOOKUP1_QUADRATIC) 
  {
    while (ABS(x - 0.5) < dx)
      x = rng_uniform_positive(rng);
  }

  if (x <= 0.5)
    assert_approx_equal(lookup1_value(table, x), x, 2.0*REAL_EPSILON);
  else 
    assert_approx_equal(lookup1_value(table, x), 0.5 - (x - 0.5), 2.0*REAL_EPSILON);

  lookup1_free(table);
}

static void test_v_linear(void** state)
{
  test_v(state, LOOKUP1_LINEAR);
}

static void test_v_quadratic(void** state)
{
  test_v(state, LOOKUP1_QUADRATIC);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_zero_linear),
    cmocka_unit_test(test_zero_quadratic),
    cmocka_unit_test(test_line_linear),
    cmocka_unit_test(test_line_quadratic),
    cmocka_unit_test(test_v_linear),
    cmocka_unit_test(test_v_quadratic)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
