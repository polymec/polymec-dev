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

#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmockery.h"
#include "core/polymec.h"
#include "core/constant_st_func.h"
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
  double o = 1.0;
  sp_func_t* one = constant_sp_func_new(1, &o);
  point_t x0 = {0.0, 0.0, 0.0};
  vector_t e3 = {.x = 0.0, .y = 0.0, .z = 1.0};
  for (int degree = 1; degree <= 5; ++degree)
  {
    sphere_integrator_t* I = sphere_integrator_new(degree);

    static const double radii[] = {1.0, 2.0, 3.0};
    for (int r = 0; r < 3; ++r)
    {
      double radius = radii[r];

      // This should give zero.
      double area;
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
static void r(void* context, point_t* x, double* result)
{
  *result = sqrt(x->x*x->x + x->y*x->y + x->z*x->z);
}

static void r2(void* context, point_t* x, double* result)
{
  *result = x->x*x->x + x->y*x->y + x->z*x->z;
}

static void r3(void* context, point_t* x, double* result)
{
  *result = pow(x->x*x->x + x->y*x->y + x->z*x->z, 1.5);
}

static void r4(void* context, point_t* x, double* result)
{
  *result = pow(x->x*x->x + x->y*x->y + x->z*x->z, 2.0);
}

static void test_cap_radial_function(void** state, 
                                     sp_func_t* f, 
                                     int f_degree,
                                     double radius, 
                                     double gamma,
                                     double answer)
{
  point_t x0 = {0.0, 0.0, 0.0};

  // Try a random orientation.
  vector_t e3;
  vector_randomize(&e3, rand, 1.0);

  sphere_integrator_t* I = sphere_integrator_new(f_degree);

  double result;
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
                                              double radius, 
                                              double gamma,
                                              double t, 
                                              double answer,
                                              double tolerance)
{
  point_t x0 = {0.0, 0.0, 0.0};

  // Try a random orientation.
  vector_t e3;
  vector_randomize(&e3, rand, 1.0);

  sphere_integrator_t* I = sphere_integrator_new(f_degree);

  double result;
  sphere_integrator_cap_at_time(I, &x0, radius, f, &e3, gamma, t, &result);
  assert_true(fabs(result - answer) < tolerance);

  sphere_integrator_free(I);
}

static void linear_growth(void* context, point_t* x, double t, double* result)
{
  *result = t * sqrt(x->x*x->x + x->y*x->y + x->z*x->z);
}

static void quadratic_growth(void* context, point_t* x, double t, double* result)
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
      double t = 1.0*i;

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
      double t = 1.0*i;

      // Integrate over the entire sphere.
      test_cap_time_dep_radial_function(state, f, degree, 2.0, M_PI, t, t*t*64.0*M_PI, 1e-10);

      // Integrate over half.
      test_cap_time_dep_radial_function(state, f, degree, 2.0, 0.5*M_PI, t, t*t*32.0*M_PI, 1e-10);
    }
  }

  f = NULL;
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
    unit_test(test_cap_quadratic_growth)
//    unit_test(test_cap_at_time)
  };
  return run_tests(tests);
}
