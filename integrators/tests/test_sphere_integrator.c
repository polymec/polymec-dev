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
#include "core/polymec.h"
#include "geometry/plane.h"
#include "integrators/sphere_integrator.h"

void test_ctor(void** state)
{
  for (int degree = 1; degree <= 5; ++degree)
  {
    sphere_integrator_t* I = sphere_integrator_new(degree);
    sphere_integrator_free(I);
  }
}

void test_cap_area(void** state)
{
  real_t o = 1.0;
  sp_func_t* one = constant_sp_func_new(1, &o);
  point_t x0 = {0.0, 0.0, 0.0};
  vector_t e3 = {.x = 0.0, .y = 0.0, .z = 1.0};
  for (int degree = 1; degree <= 5; ++degree)
  {
    sphere_integrator_t* I = sphere_integrator_new(degree);

    static const real_t radii[] = {1.0, 2.0, 3.0};
    for (int r = 0; r < 3; ++r)
    {
      real_t radius = radii[r];

      // This should give zero.
      real_t area;
      sphere_integrator_cap(I, &x0, radius, one, &e3, 0.0, &area);
      assert_true(fabs(area) < 1e-12);

      // This should give the area of the entire sphere.
      sphere_integrator_cap(I, &x0, radius, one, &e3, M_PI, &area);
      assert_true(fabs(area - 4.0*M_PI*radius*radius) < 1e-12);

      // This should give half the area of the sphere.
      sphere_integrator_cap(I, &x0, radius, one, &e3, 0.5*M_PI, &area);
      assert_true(fabs(area - 2.0*M_PI*radius*radius) < 1e-12);
    }

    sphere_integrator_free(I);
  }
  one = NULL;
}

// Radial functions.
static void r(void* context, point_t* x, real_t* result)
{
  *result = sqrt(x->x*x->x + x->y*x->y + x->z*x->z);
}

static void r2(void* context, point_t* x, real_t* result)
{
  *result = x->x*x->x + x->y*x->y + x->z*x->z;
}

static void r3(void* context, point_t* x, real_t* result)
{
  *result = pow(x->x*x->x + x->y*x->y + x->z*x->z, 1.5);
}

static void r4(void* context, point_t* x, real_t* result)
{
  *result = pow(x->x*x->x + x->y*x->y + x->z*x->z, 2.0);
}

static void test_cap_radial_function(void** state, 
                                     sp_func_t* f, 
                                     int f_degree,
                                     real_t radius, 
                                     real_t gamma,
                                     real_t answer)
{
  point_t x0 = {0.0, 0.0, 0.0};

  // Try a random orientation.
  vector_t e3;
  vector_randomize(&e3, rand, 1.0);

  sphere_integrator_t* I = sphere_integrator_new(f_degree);

  real_t result;
  sphere_integrator_cap(I, &x0, radius, f, &e3, gamma, &result);
  assert_true(fabs(result - answer) < 1e-12);

  sphere_integrator_free(I);
}

void test_cap_r(void** state)
{
  sp_func_t* f = sp_func_from_func("r", r, SP_INHOMOGENEOUS, 1);

  for (int degree = 1; degree <= 4; ++degree)
  {
    // Integrate over the entire sphere.
    test_cap_radial_function(state, f, degree, 2.0, M_PI, 32.0*M_PI);

    // Integrate over half.
    test_cap_radial_function(state, f, degree, 2.0, 0.5*M_PI, 16.0*M_PI);
  }

  f = NULL;
}

void test_cap_r2(void** state)
{
  sp_func_t* f = sp_func_from_func("r", r2, SP_INHOMOGENEOUS, 1);

  for (int degree = 2; degree <= 4; ++degree)
  {
    // Integrate over the entire sphere.
    test_cap_radial_function(state, f, degree, 2.0, M_PI, 64.0*M_PI);

    // Integrate over half.
    test_cap_radial_function(state, f, degree, 2.0, 0.5*M_PI, 32.0*M_PI);
  }

  f = NULL;
}

void test_cap_r3(void** state)
{
  sp_func_t* f = sp_func_from_func("r", r3, SP_INHOMOGENEOUS, 1);

  for (int degree = 3; degree <= 4; ++degree)
  {
    // Integrate over the entire sphere.
    test_cap_radial_function(state, f, degree, 2.0, M_PI, 128.0*M_PI);

    // Integrate over half.
    test_cap_radial_function(state, f, degree, 2.0, 0.5*M_PI, 64.0*M_PI);
  }

  f = NULL;
}

void test_cap_r4(void** state)
{
  sp_func_t* f = sp_func_from_func("r", r4, SP_INHOMOGENEOUS, 1);

  int degree = 4;
  // Integrate over the entire sphere.

  test_cap_radial_function(state, f, degree, 2.0, M_PI, 256.0*M_PI);

  // Integrate over half.
  test_cap_radial_function(state, f, degree, 2.0, 0.5*M_PI, 128.0*M_PI);

  f = NULL;
}

static void test_cap_time_dep_radial_function(void** state, 
                                              st_func_t* f, 
                                              int f_degree,
                                              real_t radius, 
                                              real_t gamma,
                                              real_t t, 
                                              real_t answer,
                                              real_t tolerance)
{
  point_t x0 = {0.0, 0.0, 0.0};

  // Try a random orientation.
  vector_t e3;
  vector_randomize(&e3, rand, 1.0);

  sphere_integrator_t* I = sphere_integrator_new(f_degree);

  real_t result;
  sphere_integrator_cap_at_time(I, &x0, radius, f, &e3, gamma, t, &result);
  assert_true(fabs(result - answer) < tolerance);

  sphere_integrator_free(I);
}

static void linear_growth(void* context, point_t* x, real_t t, real_t* result)
{
  *result = t * sqrt(x->x*x->x + x->y*x->y + x->z*x->z);
}

static void quadratic_growth(void* context, point_t* x, real_t t, real_t* result)
{
  *result = t * t * (x->x*x->x + x->y*x->y + x->z*x->z);
}

void test_cap_linear_growth(void** state)
{
  st_func_t* f = st_func_from_func("linear_growth", linear_growth, ST_INHOMOGENEOUS, ST_NONCONSTANT, 1);

  for (int degree = 1; degree <= 4; ++degree)
  {
    for (int i = 1; i <= 10; ++i)
    {
      real_t t = 1.0*i;

      // Integrate over the entire sphere.
      test_cap_time_dep_radial_function(state, f, degree, 2.0, M_PI, t, t*32.0*M_PI, 1e-12);

      // Integrate over half.
      test_cap_time_dep_radial_function(state, f, degree, 2.0, 0.5*M_PI, t, t*16.0*M_PI, 1e-12);
    }
  }

  f = NULL;
}

void test_cap_quadratic_growth(void** state)
{
  st_func_t* f = st_func_from_func("quadratic_growth", quadratic_growth, ST_INHOMOGENEOUS, ST_NONCONSTANT, 1);

  for (int degree = 2; degree <= 4; ++degree)
  {
    for (int i = 1; i <= 10; ++i)
    {
      real_t t = 1.0*i;

      // Integrate over the entire sphere.
      test_cap_time_dep_radial_function(state, f, degree, 2.0, M_PI, t, t*t*64.0*M_PI, 1e-10);

      // Integrate over half.
      test_cap_time_dep_radial_function(state, f, degree, 2.0, 0.5*M_PI, t, t*t*32.0*M_PI, 1e-10);
    }
  }

  f = NULL;
}

void test_boundary_surface_weights_for_degree_and_boundary(void** state, point_t* x0, real_t R, int degree, sp_func_t* boundary)
{
  sphere_integrator_t* I = sphere_integrator_new(degree);
  int num_weights = sphere_integrator_num_cap_points(I);
  real_t weights[num_weights];

  sphere_integrator_compute_boundary_surface_weights(I, x0, R, boundary, weights);
for (int i = 0; i < num_weights; ++i)
printf("%g ", weights[i]);
printf("\n");
}

void test_boundary_surface_weights_for_boundary(void** state, point_t* x0, real_t R, sp_func_t* boundary)
{
  for (int p = 4; p < 5; ++p)
    test_boundary_surface_weights_for_degree_and_boundary(state, x0, R, p, boundary);
}

void test_compute_boundary_surface_weights(void** state)
{
  vector_t np = {1.0, 0.0, 0.0};
  point_t xp = {0.5, 0.0, 0.0};
  point_t x0 = {0.0, 0.0, 0.0};
  sp_func_t* plane = plane_new(&np, &xp);
  real_t R = 1.0;
  test_boundary_surface_weights_for_boundary(state, &x0, R, plane);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_ctor),
    unit_test(test_cap_area),
    unit_test(test_cap_r),
    unit_test(test_cap_r2),
    unit_test(test_cap_r3),
    unit_test(test_cap_r4),
    unit_test(test_cap_linear_growth),
    unit_test(test_cap_quadratic_growth),
    unit_test(test_compute_boundary_surface_weights)
//    unit_test(test_cap_at_time)
  };
  return run_tests(tests);
}
