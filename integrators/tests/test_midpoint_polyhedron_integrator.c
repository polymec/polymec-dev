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
#include "integrators/polyhedron_integrator.h"

void test_ctor(void** state)
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

static void set_cubic_domain(polyhedron_integrator_t* I, real_t s)
{
  // We define a cube as the volume bounded by 6 faces, each defined 
  // by nodes traversed in counterclockwise order.
  point_t nodes[24] = {// -x
                       {.x = 0.0, .y = 0.0, .z = 0.0},
                       {.x = 0.0, .y = 0.0, .z = s  },
                       {.x = 0.0, .y = s,   .z = s  },
                       {.x = 0.0, .y = s,   .z = 0.0},

                       // +x
                       {.x = s,   .y = s,   .z = 0.0},
                       {.x = s,   .y = s,   .z = s  },
                       {.x = s,   .y = 0.0, .z = s  },
                       {.x = s,   .y = 0.0, .z = 0.0},

                       // -y
                       {.x = 0.0, .y = 0.0, .z = 0.0},
                       {.x = s,   .y = 0.0, .z = 0.0},
                       {.x = s,   .y = 0.0, .z = s  },
                       {.x = 0.0, .y = 0.0, .z = s  },

                       // +y
                       {.x = s,   .y = s,   .z = 0.0},
                       {.x = 0.0, .y = s,   .z = 0.0},
                       {.x = 0.0, .y = s,   .z = s  },
                       {.x = s,   .y = s,   .z = s  },

                       // -z
                       {.x = 0.0, .y = 0.0, .z = 0.0},
                       {.x = 0.0, .y = s,   .z = 0.0},
                       {.x = s,   .y = s,   .z = 0.0},
                       {.x = s,   .y = 0.0, .z = 0.0},
                       
                       // +z
                       {.x = 0.0, .y = s,   .z = s},
                       {.x = 0.0, .y = 0.0, .z = s},
                       {.x = s,   .y = 0.0, .z = s},
                       {.x = s,   .y = s,   .z = s}};
                       
  int node_offsets[7] = {0, 4, 8, 12, 16, 20, 24};
  polyhedron_integrator_set_domain(I, nodes, node_offsets, 6);
}

static real_t cube_volume_integral(polyhedron_integrator_t* I, 
                                   real_t s, 
                                   scalar_integrand_func f)
{
  set_cubic_domain(I, s);

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
  set_cubic_domain(I, s);

  real_t integral = 0.0;
  int pos = 0;
  point_t xq;
  vector_t nq;
  real_t wq;
  while (polyhedron_integrator_next_surface_point(I, &pos, &xq, &nq, &wq))
  {
    vector_t v;
    f(&xq, &v);
    integral += vector_dot(&v, &nq);
  }
  return integral;
}

void test_cube_volume_integrals(void** state)
{
  polyhedron_integrator_t* I = midpoint_polyhedron_integrator_new();
  real_t I1_1 = cube_volume_integral(I, 1.0, f1);
  assert_true(fabs(I1_1 - 0.5) < 1e-14);
  real_t I1_2 = cube_volume_integral(I, 0.5, f1);
  assert_true(fabs(I1_2 - 0.125) < 1e-14);
  real_t I2_1 = cube_volume_integral(I, 1.0, f2);
  assert_true(fabs(I2_1 - 0.75) < 1e-14);
  real_t I2_2 = cube_volume_integral(I, 0.5, f2);
  assert_true(fabs(I2_2 - 0.15625) < 1e-14);
  polyhedron_integrator_free(I);
}

void test_cube_surface_integrals(void** state)
{
  polyhedron_integrator_t* I = midpoint_polyhedron_integrator_new();
  real_t I1_1 = cube_surface_integral(I, 1.0, g1);
  assert_true(fabs(I1_1 - 1.0) < 1e-14);
  real_t I1_2 = cube_surface_integral(I, 0.5, g1);
  assert_true(fabs(I1_2 - 0.5) < 1e-14);
  real_t I2_1 = cube_surface_integral(I, 1.0, g2);
  assert_true(fabs(I2_1 - 2.0) < 1e-14);
  real_t I2_2 = cube_surface_integral(I, 0.5, g2);
  assert_true(fabs(I2_2 - 0.75) < 1e-14);
  polyhedron_integrator_free(I);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_ctor),
    unit_test(test_cube_volume_integrals),
    unit_test(test_cube_surface_integrals)
  };
  return run_tests(tests);
}

