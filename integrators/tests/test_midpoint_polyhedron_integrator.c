// Copyright (c) 2012-2016, Jeffrey N. Johnson
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
#include "core/polymec.h"
#include "geometry/create_uniform_mesh.h"
#include "integrators/polyhedron_integrator.h"

static void test_ctor(void** state)
{
  polyhedron_integrator_t* I = midpoint_polyhedron_integrator_new();
  polyhedron_integrator_free(I);
}

// Functions for integrating over cubes.
typedef real_t (*scalar_integrand_func)(point_t* x);
typedef void (*vector_integrand_func)(point_t* x, vector_t* v);

// f1(x,y,z) = x
static real_t f1(point_t* x)
{
  return x->x;
}

// f2(x,y,z) = x**2 + y
static real_t f2(point_t* x)
{
  return x->x*x->x + x->y;
}

// g1(x,y,z) = (x, 0, 0)
static void g1(point_t* x, vector_t* v)
{
  v->x = x->x;
  v->y = v->z = 0.0;
}

// g2(x,y,z) = (x**2, y)
static void g2(point_t* x, vector_t* v)
{
  v->x = x->x*x->x;
  v->y = x->y;
  v->z = 0.0;
}

static real_t cube_volume_integral(polyhedron_integrator_t* I, 
                                   real_t s, 
                                   scalar_integrand_func f)
{
  bbox_t bbox = {.x1 = 0.0, .x2 = s, .y1 = 0.0, .y2 = s, .z1 = 0.0, .z2 = s};
  mesh_t* mesh = create_uniform_mesh(MPI_COMM_SELF, 1, 1, 1, &bbox);
  polyhedron_integrator_set_domain(I, mesh, 0);
  mesh_free(mesh);

  real_t integral = 0.0;
  int pos = 0;
  point_t xq;
  real_t wq;
  while (polyhedron_integrator_next_volume_point(I, &pos, &xq, &wq))
    integral += wq * f(&xq);
  return integral;
}

static real_t cube_surface_integral(polyhedron_integrator_t* I, 
                                    real_t s, 
                                    vector_integrand_func f)
{
  bbox_t bbox = {.x1 = 0.0, .x2 = s, .y1 = 0.0, .y2 = s, .z1 = 0.0, .z2 = s};
  mesh_t* mesh = create_uniform_mesh(MPI_COMM_SELF, 1, 1, 1, &bbox);
  polyhedron_integrator_set_domain(I, mesh, 0);
  mesh_free(mesh);

  real_t integral = 0.0;
  int pos;
  point_t xq;
  vector_t nq;
  real_t wq;
  for (int face = 0; face < 6; ++face)
  {
    pos = 0;
    while (polyhedron_integrator_next_surface_point(I, face, &pos, &xq, &nq, &wq))
    {
      vector_t v;
      f(&xq, &v);
      integral += wq * vector_dot(&v, &nq);
    }
  }
  return integral;
}

static void test_cube_volume_integrals(void** state)
{
  polyhedron_integrator_t* I = midpoint_polyhedron_integrator_new();
  real_t I1_1 = cube_volume_integral(I, 1.0, f1);
  assert_true(fabs(I1_1 - 0.5) < 1e-14);
  real_t I1_2 = cube_volume_integral(I, 0.5, f1);
  assert_true(fabs(I1_2 - 0.03125) < 1e-14);
  real_t I2_1 = cube_volume_integral(I, 1.0, f2);
  assert_true(fabs(I2_1 - 0.75) < 1e-14);
  real_t I2_2 = cube_volume_integral(I, 0.5, f2);
  assert_true(fabs(I2_2 - 0.0390625) < 1e-14);
  polyhedron_integrator_free(I);
}

static void test_cube_surface_integrals(void** state)
{
  polyhedron_integrator_t* I = midpoint_polyhedron_integrator_new();
  real_t I1_1 = cube_surface_integral(I, 1.0, g1);
  assert_true(fabs(I1_1 - 1.0) < 1e-14);
  real_t I1_2 = cube_surface_integral(I, 0.5, g1);
  assert_true(fabs(I1_2 - 0.125) < 1e-14);
  real_t I2_1 = cube_surface_integral(I, 1.0, g2);
  assert_true(fabs(I2_1 - 2.0) < 1e-14);
  real_t I2_2 = cube_surface_integral(I, 0.5, g2);
  assert_true(fabs(I2_2 - 0.1875) < 1e-14);
  polyhedron_integrator_free(I);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_ctor),
    cmocka_unit_test(test_cube_volume_integrals),
    cmocka_unit_test(test_cube_surface_integrals)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}

