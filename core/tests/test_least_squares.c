// Copyright 2012-2013 Jeffrey Johnson.
// 
// This file is part of Polymec, and is licensed under the Apache License, 
// Version 2.0 (the "License"); you may not use this file except in 
// compliance with the License. You may may find the text of the license in 
// the LICENSE file at the top-level source directory, or obtain a copy of 
// it at
// 
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

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

void test_poly_shape_functions(void** state, int p, point_t* x0, point_t* points, int num_points, double* coeffs, bool weighted)
{
  polynomial_t* poly = polynomial_new(p, coeffs, x0);

  // Create scatter data.
  double data[num_points];
  for (int i = 0; i < num_points; ++i)
    data[i] = polynomial_value(poly, &points[i]);

  // Compute shape functions for the given data.
  ls_weight_func_t* W = NULL;
  if (weighted)
    W = simple_w_new(2, 1e-4);
  poly_ls_shape_t* N = poly_ls_shape_new(p, W, false);
  poly_ls_shape_set_domain(N, x0, points, num_points);
  double Nk[num_points];

  // Make sure the shape functions interpolate the data.
  for (int i = 0; i < num_points; ++i)
  {
    double value = 0.0;
    poly_ls_shape_compute(N, &points[i], Nk);
    for (int k = 0; k < num_points; ++k)
      value += Nk[k]*data[k];
    assert_true(fabs(value - data[i]) < 1e-12);
  }

  // Now make sure that the fit matches the polynomial at another point.
  point_t point;
  generate_random_points(1, &point);
  poly_ls_shape_compute(N, &point, Nk);
  double phi_fit = 0.0;
  for (int k = 0; k < num_points; ++k)
    phi_fit += Nk[k] * data[k];
  double phi = polynomial_value(poly, &point);
//  printf("%g %g %g\n", phi_fit, phi, fabs(phi_fit - phi));
  assert_true(fabs(phi_fit - phi) < 1e-12);

  // Clean up.
  N = NULL;
  poly = NULL;
}

void test_poly_shape_function_gradients(void** state, int p, point_t* x0, point_t* points, int num_points, double* coeffs, bool weighted)
{
  polynomial_t* poly = polynomial_new(p, coeffs, x0);

  // Create scatter data.
  double data[num_points];
  vector_t data_grads[num_points];
  memset(data, 0, sizeof(double)*num_points);
  memset(data_grads, 0, sizeof(vector_t)*num_points);
  for (int i = 0; i < num_points; ++i)
  {
    data[i] = polynomial_value(poly, &points[i]);
    data_grads[i].x = polynomial_deriv(poly, 1, 0, 0, &points[i]);
    data_grads[i].y = polynomial_deriv(poly, 0, 1, 0, &points[i]);
    data_grads[i].z = polynomial_deriv(poly, 0, 0, 1, &points[i]);
  }

  // Compute shape functions for the given data.
  ls_weight_func_t* W = NULL;
  if (weighted)
    W = simple_w_new(2, 1e-4);
  poly_ls_shape_t* N = poly_ls_shape_new(p, W, true);
  poly_ls_shape_set_domain(N, x0, points, num_points);

  // Make sure the shape functions interpolate the data and their gradients.
  double Nk[num_points];
  vector_t gradNk[num_points];
  for (int i = 0; i < num_points; ++i)
  {
    double value = 0.0;
    vector_t gradient = {.x = 0.0, .y = 0.0, .z = 0.0};
    poly_ls_shape_compute_gradients(N, &points[i], Nk, gradNk);
    for (int k = 0; k < num_points; ++k)
    {
      value += Nk[k]*data[k];
      gradient.x += gradNk[k].x*data[k];
      gradient.y += gradNk[k].y*data[k];
      gradient.z += gradNk[k].z*data[k];
    }
    assert_true(fabs(value - data[i]) < 2e-14);
//    printf("%g %g %g\n", gradient.x, data_grads[i].x, fabs(gradient.x - data_grads[i].x));
    assert_true(fabs(gradient.x - data_grads[i].x) < 1e-4);
    assert_true(fabs(gradient.y - data_grads[i].y) < 1e-4);
    assert_true(fabs(gradient.z - data_grads[i].z) < 1e-4);
  }

  // Now make sure that the fit matches the polynomial and its 
  // derivative at another point.
  point_t point;
  generate_random_points(1, &point);
  poly_ls_shape_compute(N, &point, Nk);
  poly_ls_shape_compute_gradients(N, &point, Nk, gradNk);
  double phi_fit = 0.0;
  vector_t grad_phi_fit = {.x = 0.0, .y = 0.0, .z = 0.0};
  for (int k = 0; k < num_points; ++k)
  {
    phi_fit += Nk[k] * data[k];
    grad_phi_fit.x += gradNk[k].x * data[k];
    grad_phi_fit.y += gradNk[k].y * data[k];
    grad_phi_fit.z += gradNk[k].z * data[k];
  }
  double phi = polynomial_value(poly, &point);
  vector_t grad_phi;
  grad_phi.x = polynomial_deriv(poly, 1, 0, 0, &point);
  grad_phi.y = polynomial_deriv(poly, 0, 1, 0, &point);
  grad_phi.z = polynomial_deriv(poly, 0, 0, 1, &point);
//  printf("%g %g %g\n", phi_fit, phi, fabs(phi_fit - phi));
  assert_true(fabs(phi_fit - phi) < 5e-14);
//printf("%g %g %g\n", grad_phi_fit.x, grad_phi.x, fabs(grad_phi_fit.x - grad_phi.x));
  assert_true(fabs(grad_phi_fit.x - grad_phi.x) < 1e-5);
  assert_true(fabs(grad_phi_fit.y - grad_phi.y) < 1e-5);
  assert_true(fabs(grad_phi_fit.z - grad_phi.z) < 1e-5);

  // Clean up.
  N = NULL;
  poly = NULL;
}

void test_poly_shape_function_constraints(void** state, int p, point_t* x0, point_t* points, int num_points, int num_ghosts, double* coeffs, bool weighted)
{
  ASSERT(p > 0);
  ASSERT(num_ghosts <= num_points/2);

  polynomial_t* poly = polynomial_new(p, coeffs, x0);

  // Create scatter data.
  double data[num_points];
  vector_t data_grads[num_points];
  memset(data, 0, sizeof(double)*num_points);
  memset(data_grads, 0, sizeof(vector_t)*num_points);
  for (int i = 0; i < num_points; ++i)
  {
    data[i] = polynomial_value(poly, &points[i]);
    data_grads[i].x = polynomial_deriv(poly, 1, 0, 0, &points[i]);
    data_grads[i].y = polynomial_deriv(poly, 0, 1, 0, &points[i]);
    data_grads[i].z = polynomial_deriv(poly, 0, 0, 1, &points[i]);
  }

  // Compute shape functions for the given data.
  ls_weight_func_t* W = NULL;
  if (weighted)
    W = simple_w_new(2, 1e-4);
  poly_ls_shape_t* N = poly_ls_shape_new(p, W, true);
  poly_ls_shape_set_domain(N, x0, points, num_points);

  // Constraints: make the constraint points match their corresponding 
  // data points by enforcing the constraint:
  // phi + 2*dphi/dx + 3*dphi/dy + 4*dhi/dz = (same thing)
  int ghost_indices[num_ghosts];
  double a[num_ghosts], b[num_ghosts], c[num_ghosts], d[num_ghosts], e[num_ghosts];
  for (int i = 0; i < num_ghosts; ++i)
  {
    int j = num_points - num_ghosts + i;
    ghost_indices[i] = j;
    a[i] = 1.0, b[i] = 2.0, c[i] = 3.0, d[i] = 4.0; 
    e[i] = data[j] + 2.0*data_grads[j].x + 3.0*data_grads[j].y + 4.0*data_grads[j].z;
  }

  // Compute the last (num_ghosts) points using the first 
  // (num_points - num_ghosts) points and constraints. For 
  // this, we construct the affine transformation that maps the 
  // total set of solution values to the constrained values.
  double A[num_ghosts*num_points], B[num_ghosts];
  point_t xc[num_ghosts];
  for (int c = 0; c < num_ghosts; ++c)
  {
    xc[c].x = points[ghost_indices[c]].x;
    xc[c].y = points[ghost_indices[c]].y;
    xc[c].z = points[ghost_indices[c]].z;
  }
  poly_ls_shape_compute_ghost_transform(N, ghost_indices, num_ghosts, xc,
                                        a, b, c, d, e, A, B);

  // Now apply the affine transformation to the data and make sure we 
  // reproduce the data points. Remember that A is stored in column-
  // major order!
  double constrained_data[num_ghosts];
  for (int i = 0; i < num_ghosts; ++i)
  {
    constrained_data[i] = B[i];
    for (int j = 0; j < num_points; ++j)
      constrained_data[i] += A[num_ghosts*j+i] * data[j];
//printf("%g %g %g\n", constrained_data[i], data[ghost_indices[i]], fabs(constrained_data[i] - data[ghost_indices[i]]));
    assert_true(fabs(constrained_data[i] - data[ghost_indices[i]]) < 8e-10);
  }

  // Clean up.
  N = NULL;
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

void test_p0_shape_funcs(void** state)
{
  static double coeffs[] = {1.0};
  point_t x0, points[4];
  generate_random_points(4, points);
  average_points(points, 4, &x0);
  test_poly_shape_functions(state, 0, &x0, points, 4, coeffs, false);
}

void test_weighted_p0_shape_funcs(void** state)
{
  static double coeffs[] = {1.0};
  point_t x0, points[4];
  generate_random_points(4, points);
  average_points(points, 4, &x0);
  test_poly_shape_functions(state, 0, &x0, points, 4, coeffs, true);
}

void test_p0_shape_func_gradients(void** state)
{
  static double coeffs[] = {1.0};
  point_t x0, points[4];
  generate_random_points(4, points);
  average_points(points, 4, &x0);
  test_poly_shape_function_gradients(state, 0, &x0, points, 4, coeffs, false);
}

void test_weighted_p0_shape_func_gradients(void** state)
{
  static double coeffs[] = {1.0};
  point_t x0, points[4];
  generate_random_points(4, points);
  average_points(points, 4, &x0);
  test_poly_shape_function_gradients(state, 0, &x0, points, 4, coeffs, true);
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

void test_p1_shape_funcs(void** state)
{
  static double coeffs[] = {1.0, 2.0, 3.0, 4.0};
  point_t x0, points[8];
  generate_random_points(8, points);
  average_points(points, 8, &x0);
  test_poly_shape_functions(state, 1, &x0, points, 8, coeffs, false);
}

void test_weighted_p1_shape_funcs(void** state)
{
  static double coeffs[] = {1.0, 2.0, 3.0, 4.0};
  point_t x0, points[8];
  generate_random_points(8, points);
  average_points(points, 8, &x0);
  test_poly_shape_functions(state, 1, &x0, points, 8, coeffs, true);
}

void test_p1_shape_func_gradients(void** state)
{
  static double coeffs[] = {1.0, 2.0, 3.0, 4.0};
  point_t x0, points[8];
  generate_random_points(8, points);
  average_points(points, 8, &x0);
  test_poly_shape_function_gradients(state, 1, &x0, points, 8, coeffs, false);
}

void test_weighted_p1_shape_func_gradients(void** state)
{
  static double coeffs[] = {1.0, 2.0, 3.0, 4.0};
  point_t x0, points[8];
  generate_random_points(8, points);
  average_points(points, 8, &x0);
  test_poly_shape_function_gradients(state, 1, &x0, points, 8, coeffs, true);
}

void test_p1_shape_func_constraints(void** state)
{
  static double coeffs[] = {1.0, 2.0, 3.0, 4.0};
  point_t x0, points[8];
  generate_random_points(8, points);
  average_points(points, 8, &x0);
  test_poly_shape_function_constraints(state, 1, &x0, points, 8, 2, coeffs, false);
}

void test_weighted_p1_shape_func_constraints(void** state)
{
  static double coeffs[] = {1.0, 2.0, 3.0, 4.0};
  point_t x0, points[8];
  generate_random_points(8, points);
  average_points(points, 8, &x0);
  test_poly_shape_function_constraints(state, 1, &x0, points, 8, 2, coeffs, true);
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

void test_p2_shape_funcs(void** state)
{
  static double coeffs[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
  point_t x0, points[16];
  generate_random_points(16, points);
  average_points(points, 16, &x0);
  test_poly_shape_functions(state, 2, &x0, points, 16, coeffs, false);
}

void test_weighted_p2_shape_funcs(void** state)
{
  static double coeffs[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
  point_t x0, points[16];
  generate_random_points(16, points);
  average_points(points, 16, &x0);
  test_poly_shape_functions(state, 2, &x0, points, 16, coeffs, true);
}

void test_p2_shape_func_gradients(void** state)
{
  static double coeffs[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
  point_t x0, points[16];
  generate_random_points(16, points);
  average_points(points, 16, &x0);
  test_poly_shape_function_gradients(state, 2, &x0, points, 16, coeffs, false);
}

void test_weighted_p2_shape_func_gradients(void** state)
{
  static double coeffs[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
  point_t x0, points[16];
  generate_random_points(16, points);
  average_points(points, 16, &x0);
  test_poly_shape_function_gradients(state, 2, &x0, points, 16, coeffs, true);
}

void test_p2_shape_func_constraints(void** state)
{
  static double coeffs[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
  point_t x0, points[16];
  generate_random_points(16, points);
  average_points(points, 16, &x0);
  test_poly_shape_function_constraints(state, 2, &x0, points, 16, 6, coeffs, false);
}

void test_weighted_p2_shape_func_constraints(void** state)
{
  static double coeffs[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
  point_t x0, points[16];
  generate_random_points(16, points);
  average_points(points, 16, &x0);
  test_poly_shape_function_constraints(state, 2, &x0, points, 16, 6, coeffs, true);
}

void test_p3_fit(void** state)
{
  static double coeffs[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
  point_t x0, points[30];
  generate_random_points(30, points);
  test_poly_fit(state, 3, NULL, points, 30, coeffs, false);
  average_points(points, 30, &x0);
  test_poly_fit(state, 3, &x0, points, 30, coeffs, false);
}

void test_weighted_p3_fit(void** state)
{
  static double coeffs[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
  point_t x0, points[30];
  generate_random_points(30, points);
  test_poly_fit(state, 3, NULL, points, 30, coeffs, true);
  average_points(points, 30, &x0);
  test_poly_fit(state, 3, &x0, points, 30, coeffs, true);
}

void test_p3_shape_funcs(void** state)
{
  static double coeffs[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
  point_t x0, points[30];
  generate_random_points(30, points);
  average_points(points, 30, &x0);
  test_poly_shape_functions(state, 3, &x0, points, 30, coeffs, false);
}

void test_weighted_p3_shape_funcs(void** state)
{
  static double coeffs[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
  point_t x0, points[30];
  generate_random_points(30, points);
  average_points(points, 30, &x0);
  test_poly_shape_functions(state, 3, &x0, points, 30, coeffs, true);
}

void test_p3_shape_func_gradients(void** state)
{
  static double coeffs[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
  point_t x0, points[30];
  generate_random_points(30, points);
  average_points(points, 30, &x0);
  test_poly_shape_function_gradients(state, 3, &x0, points, 30, coeffs, false);
}

void test_weighted_p3_shape_func_gradients(void** state)
{
  static double coeffs[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
  point_t x0, points[30];
  generate_random_points(30, points);
  average_points(points, 30, &x0);
  test_poly_shape_function_gradients(state, 3, &x0, points, 30, coeffs, true);
}

void test_p3_shape_func_constraints(void** state)
{
  static double coeffs[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
  point_t x0, points[30];
  generate_random_points(30, points);
  average_points(points, 30, &x0);
  test_poly_shape_function_constraints(state, 3, &x0, points, 30, 12, coeffs, false);
}

void test_weighted_p3_shape_func_constraints(void** state)
{
  static double coeffs[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
  point_t x0, points[30];
  generate_random_points(30, points);
  average_points(points, 30, &x0);
  test_poly_shape_function_constraints(state, 3, &x0, points, 30, 12, coeffs, true);
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
    unit_test(test_weighted_p3_fit),
    unit_test(test_p0_shape_funcs),
    unit_test(test_p1_shape_funcs),
    unit_test(test_p2_shape_funcs),
    unit_test(test_p3_shape_funcs),
    unit_test(test_weighted_p0_shape_funcs),
    unit_test(test_weighted_p1_shape_funcs),
    unit_test(test_weighted_p2_shape_funcs),
    unit_test(test_weighted_p3_shape_funcs)
    // We don't need this stuff at the moment.
//    unit_test(test_p0_shape_func_gradients),
//    unit_test(test_p1_shape_func_gradients),
//    unit_test(test_p2_shape_func_gradients),
//    unit_test(test_p3_shape_func_gradients),
//    unit_test(test_weighted_p0_shape_func_gradients),
//    unit_test(test_weighted_p1_shape_func_gradients),
//    unit_test(test_weighted_p2_shape_func_gradients),
//    unit_test(test_weighted_p3_shape_func_gradients),
//    unit_test(test_p1_shape_func_constraints),
//    unit_test(test_p2_shape_func_constraints),
//    unit_test(test_p3_shape_func_constraints),
//    unit_test(test_weighted_p1_shape_func_constraints),
//    unit_test(test_weighted_p2_shape_func_constraints),
//    unit_test(test_weighted_p3_shape_func_constraints)
  };
  return run_tests(tests);
}
