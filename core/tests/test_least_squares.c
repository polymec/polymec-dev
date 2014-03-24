// Copyright (c) 2012-2014, Jeffrey N. Johnson
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this 
// list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice, 
// this list of conditions and the following disclaimer in the documentation 
// and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include <stdlib.h>
#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include <time.h>
#include "cmockery.h"
#include "core/least_squares.h"
#include "core/linear_algebra.h"

static void generate_random_points(int num_points, point_t* points)
{
  static bbox_t bbox = {.x1 = -1.0, .x2 = 1.0, .y1 = -1.0, .y2 = 1.0, .z1 = -1.0, .z2 = 1.0};

  for (int i = 0; i < num_points; ++i)
    point_randomize(&points[i], rand, &bbox);
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
  double A, B, h;
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

static void simple_W_eval(void* context, vector_t* y, double* W, vector_t* gradient)
{
  simple_W* wf = context;
  double D = vector_mag(y) / wf->h;
  *W = 1.0 / (pow(D, wf->A) + pow(wf->B, wf->A));
  if (D == 0.0)
  {
    gradient->x = gradient->y = gradient->z = 0.0;
  }
  else
  {
    double dDdx = y->x / D, dDdy = y->y / D, dDdz = y->z / D;
    double deriv_term = -(*W)*(*W) * 2.0 * D;
    gradient->x = deriv_term * dDdx;
    gradient->y = deriv_term * dDdy;
    gradient->z = deriv_term * dDdz;
  }
}

static ls_weight_func_t* simple_w_new(int A, double B)
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

void test_poly_fit(void** state, int p, point_t* x0, point_t* points, int num_points, double* coeffs, bool weighted)
{
  polynomial_t* poly = polynomial_new(p, coeffs, x0);

  // Create scatter data using a polynomial with the given coefficients.
  double data[num_points];
  for (int i = 0; i < num_points; ++i)
    data[i] = polynomial_value(poly, &points[i]);

  // Assemble the linear system to recover the coefficients.
  int dim = polynomial_num_terms(poly);
  double A[dim*dim], b[dim];
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
  dgetrf(&dim, &dim, A, &lda, pivot, &info);
  assert_int_equal(0, info);
  dgetrs(&trans, &dim, &one, A, &lda, pivot, b, &ldb, &info);
  assert_int_equal(0, info);

  // Make sure we've recovered the coefficients.
  for (int i = 0; i < dim; ++i)
  {
//  printf("%g %g %g\n", b[i], coeffs[i], fabs(b[i] - coeffs[i]));
// FIXME: This fails sometimes.
    assert_true(fabs(b[i] - coeffs[i]) < 1e-12);
  }

  poly = NULL;
}

void test_p0_fit(void** state)
{
  static double coeffs[] = {1.0};
  point_t x0, points[4];
  generate_random_points(4, points);
  test_poly_fit(state, 0, NULL, points, 4, coeffs, false);
  average_points(points, 4, &x0);
  test_poly_fit(state, 0, &x0, points, 4, coeffs, false);
}

void test_weighted_p0_fit(void** state)
{
  static double coeffs[] = {1.0};
  point_t x0, points[4];
  generate_random_points(4, points);
  test_poly_fit(state, 0, NULL, points, 4, coeffs, true);
  average_points(points, 4, &x0);
  test_poly_fit(state, 0, &x0, points, 4, coeffs, true);
}

void test_p1_fit(void** state)
{
  static double coeffs[] = {1.0, 2.0, 3.0, 4.0};
  point_t x0, points[8];
  generate_random_points(8, points);
  test_poly_fit(state, 1, NULL, points, 8, coeffs, false);
  average_points(points, 8, &x0);
  test_poly_fit(state, 1, &x0, points, 8, coeffs, false);
}

void test_weighted_p1_fit(void** state)
{
  static double coeffs[] = {1.0, 2.0, 3.0, 4.0};
  point_t x0, points[8];
  generate_random_points(8, points);
  test_poly_fit(state, 1, NULL, points, 8, coeffs, true);
  average_points(points, 8, &x0);
  test_poly_fit(state, 1, &x0, points, 8, coeffs, true);
}

void test_p2_fit(void** state)
{
  static double coeffs[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
  point_t x0, points[16];
  generate_random_points(16, points);
  test_poly_fit(state, 2, NULL, points, 16, coeffs, false);
  average_points(points, 16, &x0);
  test_poly_fit(state, 2, &x0, points, 16, coeffs, false);
}

void test_weighted_p2_fit(void** state)
{
  static double coeffs[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
  point_t x0, points[16];
  generate_random_points(16, points);
  test_poly_fit(state, 2, NULL, points, 16, coeffs, true);
  average_points(points, 16, &x0);
  test_poly_fit(state, 2, &x0, points, 16, coeffs, true);
}

void test_p3_fit(void** state)
{
  static double coeffs[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 
                            11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0};
  point_t x0, points[30];
  generate_random_points(30, points);
  test_poly_fit(state, 3, NULL, points, 30, coeffs, false);
  average_points(points, 30, &x0);
  test_poly_fit(state, 3, &x0, points, 30, coeffs, false);
}

void test_weighted_p3_fit(void** state)
{
  static double coeffs[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 
                            11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0};
  point_t x0, points[30];
  generate_random_points(30, points);
  test_poly_fit(state, 3, NULL, points, 30, coeffs, true);
  average_points(points, 30, &x0);
  test_poly_fit(state, 3, &x0, points, 30, coeffs, true);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);

  // Initialize the random number generator.
  srand((unsigned)time(NULL));

  const UnitTest tests[] = 
  {
    unit_test(test_p0_fit),
    unit_test(test_p1_fit),
    unit_test(test_p2_fit), 
    unit_test(test_p3_fit),
    unit_test(test_weighted_p0_fit),
    unit_test(test_weighted_p1_fit),
    unit_test(test_weighted_p2_fit), 
    unit_test(test_weighted_p3_fit)
  };
  return run_tests(tests);
}
