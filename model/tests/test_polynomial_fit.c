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

#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmockery.h"
#include "core/polynomial.h"
#include "model/polynomial_fit.h"

void test_polynomial_fit_ctor(void** state, int num_components, int p)
{
  polynomial_fit_t* fit = polynomial_fit_new(num_components, p);
  assert_int_equal(num_components, polynomial_fit_num_components(fit));
  assert_int_equal(p, polynomial_fit_degree(fit));
  assert_int_equal(0, polynomial_fit_num_equations(fit));
  polynomial_fit_free(fit);
}

void test_polynomial_fit_new(void** state)
{
  for (int num_comp = 1; num_comp < 8; ++num_comp)
  {
    for (int degree = 0; degree < 4; ++degree)
      test_polynomial_fit_ctor(state, num_comp, degree);
  }
}

void test_fit_consistency(void** state, polynomial_fit_t* fit, int component,
                          polynomial_t* p, int num_points, bbox_t* bbox)
{
  rng_t* rng = host_rng_new();
  int num_comp = polynomial_fit_num_components(fit);

  for (int i = 0; i < num_points; ++i)
  {
    point_t xrand;
    point_randomize(&xrand, rng, bbox);

    // Polynomial values.
    real_t q = polynomial_value(p, &xrand);
    real_t dqdx = polynomial_deriv_value(p, 1, 0, 0, &xrand);
    real_t dqdy = polynomial_deriv_value(p, 0, 1, 0, &xrand);
    real_t dqdz = polynomial_deriv_value(p, 0, 0, 1, &xrand);

    // Polynomial fit values.
    real_t u[num_comp], dudx[num_comp], dudy[num_comp], dudz[num_comp];
    polynomial_fit_eval(fit, &xrand, u);
    polynomial_fit_eval_deriv(fit, 1, 0, 0, &xrand, dudx);
    polynomial_fit_eval_deriv(fit, 0, 1, 0, &xrand, dudy);
    polynomial_fit_eval_deriv(fit, 0, 0, 1, &xrand, dudz);

//printf("u error: %g-%g=%g\n", u[component], q, fabs(u[component]-q));
//printf("dudx error: %g\n", fabs(dudx[component]-dqdx));
//printf("dudy error: %g\n", fabs(dudy[component]-dqdy));
//printf("dudz error: %g\n", fabs(dudz[component]-dqdz));
    assert_true(fabs(u[component] - q) < 5e-15);
    assert_true(fabs(dudx[component] - dqdx) < 5e-15);
    assert_true(fabs(dudy[component] - dqdy) < 5e-15);
    assert_true(fabs(dudz[component] - dqdz) < 5e-15);
  }
}

void test_polynomial_fit(void** state, polynomial_t** polynomials, 
                         int num_components)
{
  rng_t* rng = host_rng_new();
  int p = polynomial_degree(polynomials[0]);

  polynomial_fit_t* fit = polynomial_fit_new(num_components, p);
  point_t x[] = {{.x =  0.0, .y =  0.0, .z =  0.0},
                 {.x =  0.0, .y =  0.0, .z = -0.5},
                 {.x =  0.0, .y = -0.5, .z =  0.0},
                 {.x = -0.5, .y =  0.0, .z =  0.0},
                 {.x =  0.0, .y =  0.0, .z =  0.5},
                 {.x =  0.0, .y =  0.5, .z =  0.0},
                 {.x =  0.5, .y =  0.0, .z =  0.0}};

  // Fit centered at x0 = {0, 0, 0}.
  point_t x0 = {.x = 0.0, .y = 0.0, .z = 0.0};
  for (int c = 0; c < num_components; ++c)
    *polynomial_x0(polynomials[c]) = x0;

  for (int i = 0; i < 7; ++i)
  {
    for (int c = 0; c < num_components; ++c)
      polynomial_fit_add_scatter_datum(fit, c, polynomial_value(polynomials[c], &x[i]), &x[i]);
  }
  assert_int_equal(7*num_components, polynomial_fit_num_equations(fit));
  polynomial_fit_compute(fit);
  bbox_t bbox = {.x1 = -0.5, .y1 = -0.5, .z1 = -0.5, .x2 = 0.5, .y2 = 0.5, .z2 = 0.5};
  for (int c = 0; c < num_components; ++c)
    test_fit_consistency(state, fit, c, polynomials[c], 7, &bbox);

  // Fit centered at a random point.
  point_randomize(&x0, rng, &bbox);
  for (int c = 0; c < num_components; ++c)
    *polynomial_x0(polynomials[c]) = x0;
  polynomial_fit_reset(fit, &x0);
  bbox.x1 += x0.x, bbox.x2 += x0.x;
  bbox.y1 += x0.y, bbox.y2 += x0.y;
  bbox.z1 += x0.z, bbox.z2 += x0.z;

  assert_int_equal(0, polynomial_fit_num_equations(fit));
  for (int i = 0; i < 7; ++i)
  {
    for (int c = 0; c < num_components; ++c)
      polynomial_fit_add_scatter_datum(fit, c, polynomial_value(polynomials[c], &x[i]), &x[i]);
  }
  polynomial_fit_compute(fit);
  for (int c = 0; c < num_components; ++c)
    test_fit_consistency(state, fit, c, polynomials[c], 7, &bbox);

  polynomial_fit_free(fit);
}

void test_polynomial_fit_constant_scalar(void** state)
{
  real_t coeffs[] = {3.0};
  point_t x0 = {.x = 0.0, .y = 0.0, .z = 0.0};
  polynomial_t* p0 = polynomial_new(0, coeffs, &x0); 
  test_polynomial_fit(state, &p0, 1);
}

void test_polynomial_fit_constant_vector(void** state)
{
  real_t coeffs1[] = {1.0};
  real_t coeffs2[] = {2.0};
  real_t coeffs3[] = {3.0};
  point_t x0 = {.x = 0.0, .y = 0.0, .z = 0.0};
  polynomial_t* p1 = polynomial_new(0, coeffs1, &x0); 
  polynomial_t* p2 = polynomial_new(0, coeffs2, &x0); 
  polynomial_t* p3 = polynomial_new(0, coeffs3, &x0); 
  polynomial_t* poly[] = {p1, p2, p3};
  test_polynomial_fit(state, poly, 3);
}

void test_polynomial_fit_linear_scalar(void** state)
{
  real_t coeffs[] = {1.0, 2.0, 3.0, 4.0};
  point_t x0 = {.x = 0.0, .y = 0.0, .z = 0.0};
  polynomial_t* p1 = polynomial_new(1, coeffs, &x0); 
  test_polynomial_fit(state, &p1, 1);
}

void test_polynomial_fit_linear_vector(void** state)
{
  real_t coeffs[] = {1.0, 2.0, 3.0, 4.0};
  point_t x0 = {.x = 0.0, .y = 0.0, .z = 0.0};
  polynomial_t* p1 = polynomial_new(1, coeffs, &x0); 
  polynomial_t* p2 = polynomial_new(1, coeffs, &x0); 
  polynomial_t* p3 = polynomial_new(1, coeffs, &x0); 
  polynomial_t* poly[] = {p1, p2, p3};
  test_polynomial_fit(state, poly, 3);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_polynomial_fit_new),
    unit_test(test_polynomial_fit_constant_scalar),
    unit_test(test_polynomial_fit_constant_vector),
    unit_test(test_polynomial_fit_linear_scalar),
    unit_test(test_polynomial_fit_linear_vector)
  };
  return run_tests(tests);
}
