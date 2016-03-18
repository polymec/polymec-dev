// Copyright (c) 2012-2016, Jeffrey N. Johnson
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
#include <time.h>
#include "cmocka.h"
#include "core/least_squares.h"
#include "core/linear_algebra.h"
#include "core/polynomial.h"

static void generate_random_points(int num_points, point_t* points)
{
  static bbox_t bbox = {.x1 = -1.0, .x2 = 1.0, .y1 = -1.0, .y2 = 1.0, .z1 = -1.0, .z2 = 1.0};

  rng_t* rng = host_rng_new();
  for (int i = 0; i < num_points; ++i)
    point_randomize(&points[i], rng, &bbox);
}

static void average_points(point_t* points, int num_points, point_t* average)
{
  average->x = 0.0;
  average->y = 0.0;
  average->z = 0.0;
  for (int i = 0; i < num_points; ++i)
  {
    average->x += points[i].x;
    average->y += points[i].y;
    average->z += points[i].z;
  }
  average->x /= num_points;
  average->y /= num_points;
  average->z /= num_points;
}

// Simple weighting function: W(x, x0) = 1/((|x-x0|/h) + B)^A
typedef struct
{
  real_t A, B, h;
} simple_W;

static void simple_W_set_domain(void* context, point_t* x0, point_t* points, int num_points)
{
  simple_W* wf = context;

  // Compute the average distance between the points. This will serve as 
  // our spatial scale length, h.
  wf->h = 0.0;
  for (int n = 0; n < num_points; ++n)
    for (int l = n+1; l < num_points; ++l)
      wf->h += point_distance(&points[n], &points[l]);
  wf->h /= (num_points*(num_points+1)/2 - num_points);
}

static void simple_W_eval(void* context, vector_t* y, real_t* W, vector_t* gradient)
{
  simple_W* wf = context;
  real_t D = vector_mag(y) / wf->h;
  *W = 1.0 / (pow(D, wf->A) + pow(wf->B, wf->A));
  if (D == 0.0)
  {
    gradient->x = gradient->y = gradient->z = 0.0;
  }
  else
  {
    real_t dDdx = y->x / D, dDdy = y->y / D, dDdz = y->z / D;
    real_t deriv_term = -(*W)*(*W) * 2.0 * D;
    gradient->x = deriv_term * dDdx;
    gradient->y = deriv_term * dDdy;
    gradient->z = deriv_term * dDdz;
  }
}

static ls_weight_func_t* simple_w_new(int A, real_t B)
{
  ASSERT(A > 0);
  ASSERT(B > 0.0);

  simple_W* W_data = malloc(sizeof(simple_W));
  W_data->A = A;
  W_data->B = B;
  ls_weight_func_vtable vtable = {.set_domain = simple_W_set_domain,
                                  .eval = simple_W_eval,
                                  .dtor = free};
  return ls_weight_func_new("Simple", W_data, vtable);
}

bool test_poly_fit(void** state, int p, point_t* x0, point_t* points, int num_points, real_t* coeffs, bool weighted)
{
  polynomial_t* poly = polynomial_new(p, coeffs, x0);

  // Create scatter data using a polynomial with the given coefficients.
  real_t data[num_points];
  for (int i = 0; i < num_points; ++i)
    data[i] = polynomial_value(poly, &points[i]);

  // Assemble the linear system to recover the coefficients.
  int dim = polynomial_num_terms(poly);
  real_t A[dim*dim], b[dim];
  if (weighted)
  {
    ls_weight_func_t* W = simple_w_new(2, 1e-4);
    compute_weighted_poly_ls_system(p, W, x0, points, num_points, data, A, b);
    W = NULL;
  }
  else
    compute_poly_ls_system(p, x0, points, num_points, data, A, b);

  // Solve the system.
  int lda = dim, ldb = dim, pivot[dim], info, one = 1;
  char trans = 'N';
  rgetrf(&dim, &dim, A, &lda, pivot, &info);
  assert_int_equal(0, info);
  rgetrs(&trans, &dim, &one, A, &lda, pivot, b, &ldb, &info);
  assert_int_equal(0, info);

  // Pick a random point (x_star) near the geometric mean (x_bar) of our points. 
  point_t x_bar = {.x = 0.0, .y = 0.0, .z = 0.0};
  for (int i = 0; i < num_points; ++i)
  {
    x_bar.x += points[i].x;
    x_bar.y += points[i].y;
    x_bar.z += points[i].z;
  }
  x_bar.x /= num_points;
  x_bar.y /= num_points;
  x_bar.z /= num_points;
  real_t L = -FLT_MAX;
  for (int i = 0; i < num_points; ++i)
    L = MAX(L, MAX(ABS(points[i].x - x_bar.x), MAX(ABS(points[i].y - x_bar.y), ABS(points[i].z - x_bar.z))));
  vector_t dx;
  rng_t* rng = host_rng_new();
  vector_randomize(&dx, rng, 0.2*L);
  point_t x_star = {.x = x_bar.x + dx.x, .y = x_bar.y + dx.y, .z = x_bar.z + dx.z};

  // Compute the "actual" and fitted values at x_star, and make sure we 
  // recover the approximation to within the accuracy of the linear
  // solve (say ~5e-12).
  real_t actual_value = polynomial_value(poly, &x_star);
  polynomial_t* fitted_poly = polynomial_new(p, b, x0);
  real_t fitted_value = polynomial_value(fitted_poly, &x_star);
  real_t error = ABS(fitted_value - actual_value);
  log_debug("error for order %d fit at x = (%g, %g, %g) with %d points: %g", 
            p, x_star.x, x_star.y, x_star.z, num_points, error); 
  bool passed;
  if (weighted)
    passed = (error < 5e-12);
  else
    passed = (error < 1e-11);

  poly = NULL;
  fitted_poly = NULL;

  return passed;
}

// We want 99% reliability for randomly selected points.

void test_p0_fit(void** state)
{
  int num_failures = 0;
  for (int i = 0; i < 100; ++i) 
  {
    static real_t coeffs[] = {1.0};
    point_t x0, points[4];
    generate_random_points(4, points);
    bool passed = test_poly_fit(state, 0, NULL, points, 4, coeffs, false);
    average_points(points, 4, &x0);
    bool passed_x0 = test_poly_fit(state, 0, &x0, points, 4, coeffs, false);

    if (!passed || !passed_x0)
      ++num_failures;
  }
  assert_true(num_failures <= 1);
}

void test_weighted_p0_fit(void** state)
{
  int num_failures = 0;
  for (int i = 0; i < 100; ++i) 
  {
    static real_t coeffs[] = {1.0};
    point_t x0, points[4];
    generate_random_points(4, points);
    bool passed = test_poly_fit(state, 0, NULL, points, 4, coeffs, true);
    average_points(points, 4, &x0);
    bool passed_x0 = test_poly_fit(state, 0, &x0, points, 4, coeffs, true);

    if (!passed || !passed_x0)
      ++num_failures;
  }
  assert_true(num_failures <= 1);
}

void test_p1_fit(void** state)
{
  int num_failures = 0;
  for (int i = 0; i < 100; ++i) 
  {
    static real_t coeffs[] = {1.0, 2.0, 3.0, 4.0};
    point_t x0, points[8];
    generate_random_points(8, points);
    bool passed = test_poly_fit(state, 1, NULL, points, 8, coeffs, false);
    average_points(points, 8, &x0);
    bool passed_x0 = test_poly_fit(state, 1, &x0, points, 8, coeffs, false);

    if (!passed || !passed_x0)
      ++num_failures;
  }
  assert_true(num_failures <= 1);
}

void test_weighted_p1_fit(void** state)
{
  int num_failures = 0;
  for (int i = 0; i < 100; ++i) 
  {
    static real_t coeffs[] = {1.0, 2.0, 3.0, 4.0};
    point_t x0, points[8];
    generate_random_points(8, points);
    bool passed = test_poly_fit(state, 1, NULL, points, 8, coeffs, true);
    average_points(points, 8, &x0);
    bool passed_x0 = test_poly_fit(state, 1, &x0, points, 8, coeffs, true);

    if (!passed || !passed_x0)
      ++num_failures;
  }
  assert_true(num_failures <= 1);
}

void test_p2_fit(void** state)
{
  int num_failures = 0;
  for (int i = 0; i < 100; ++i) 
  {
    static real_t coeffs[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
    point_t x0, points[16];
    generate_random_points(16, points);
    bool passed = test_poly_fit(state, 2, NULL, points, 16, coeffs, false);
    average_points(points, 16, &x0);
    bool passed_x0 = test_poly_fit(state, 2, &x0, points, 16, coeffs, false);

    if (!passed || !passed_x0)
      ++num_failures;
  }
  assert_true(num_failures <= 1);
}

void test_weighted_p2_fit(void** state)
{
  int num_failures = 0;
  for (int i = 0; i < 100; ++i) 
  {
    static real_t coeffs[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
    point_t x0, points[16];
    generate_random_points(16, points);
    bool passed = test_poly_fit(state, 2, NULL, points, 16, coeffs, true);
    average_points(points, 16, &x0);
    bool passed_x0 = test_poly_fit(state, 2, &x0, points, 16, coeffs, true);

    if (!passed || !passed_x0)
      ++num_failures;
  }
  assert_true(num_failures <= 1);
}

void test_p3_fit(void** state)
{
  int num_failures = 0;
  for (int i = 0; i < 100; ++i) 
  {
    static real_t coeffs[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 
      11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0};
    point_t x0, points[40];
    generate_random_points(40, points);
    bool passed = test_poly_fit(state, 3, NULL, points, 40, coeffs, false);
    average_points(points, 40, &x0);
    bool passed_x0 = test_poly_fit(state, 3, &x0, points, 40, coeffs, false);

    if (!passed || !passed_x0)
      ++num_failures;
  }
  assert_true(num_failures <= 1);
}

void test_weighted_p3_fit(void** state)
{
  int num_failures = 0;
  for (int i = 0; i < 100; ++i) 
  {
    static real_t coeffs[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 
      11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0};
    point_t x0, points[50];
    generate_random_points(50, points);
    bool passed = test_poly_fit(state, 3, NULL, points, 50, coeffs, true);
    average_points(points, 50, &x0);
    bool passed_x0 = test_poly_fit(state, 3, &x0, points, 50, coeffs, true);

    if (!passed || !passed_x0)
      ++num_failures;
  }
  assert_true(num_failures <= 1);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  set_log_level(LOG_DEBUG);

  // Initialize the random number generator.
  srand((unsigned)time(NULL));

  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_p0_fit),
    cmocka_unit_test(test_p1_fit),
    cmocka_unit_test(test_p2_fit), 
    cmocka_unit_test(test_p3_fit),
    cmocka_unit_test(test_weighted_p0_fit),
    cmocka_unit_test(test_weighted_p1_fit),
    cmocka_unit_test(test_weighted_p2_fit), 
    cmocka_unit_test(test_weighted_p3_fit)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
