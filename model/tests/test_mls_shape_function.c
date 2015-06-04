// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmockery.h"
#include "core/partition_point_cloud.h"
#include "geometry/create_point_lattice.h"
#include "model/mls_shape_function.h"

// Helper function to construct lattices of points, stencils, h fields.
static void make_lattice(int nx, int ny, int nz, real_t h_over_dx,
                         point_cloud_t** domain,
                         stencil_t** neighborhoods,
                         real_t** smoothing_lengths)
{
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  *domain = create_uniform_point_lattice(MPI_COMM_SELF, nx, ny, nz, &bbox);
  exchanger_t* ex = partition_point_cloud(domain, MPI_COMM_WORLD, NULL, 1.05);
  exchanger_free(ex);
  int num_local_points = (*domain)->num_points;
  real_t dx = 1.0 / nx;

  // Set up a "radius" field to measure point extents.
  real_t R[num_local_points];
  for (int i = 0; i < num_local_points; ++i)
    R[i] = 2.0 * h_over_dx * dx;

  // Create the stencil.
  *neighborhoods = distance_based_point_stencil_new(*domain, R);

  // Now create the smoothing lengths field.
  int num_remote_points = (*domain)->num_ghosts;
  *smoothing_lengths = polymec_malloc(sizeof(real_t) * (num_local_points + num_remote_points));
  for (int i = 0; i < num_local_points; ++i)
    (*smoothing_lengths)[i] = 0.5 * R[i];
}

void test_mls_shape_function_ctor(void** state, int p)
{
  point_cloud_t* domain;
  stencil_t* neighborhoods;
  real_t* smoothing_lengths;
  make_lattice(10, 10, 10, 1.0, &domain, &neighborhoods, &smoothing_lengths);

  shape_function_kernel_t* W = simple_shape_function_kernel_new();
  shape_function_t* phi = mls_shape_function_new(p, W, domain, neighborhoods, smoothing_lengths);

  // Try setting the neighborhoods.
  for (int i = 0; i < domain->num_points; ++i)
    shape_function_set_neighborhood(phi, i);

  // Clean up.
  shape_function_free(phi);
  point_cloud_free(domain);
  stencil_free(neighborhoods);
  polymec_free(smoothing_lengths);
}

void test_mls_shape_function_ctors(void** state)
{
  for (int p = 0; p < 4; ++p)
    test_mls_shape_function_ctor(state, p);
}

void test_mls_shape_function_zero_consistency_p(void** state, int p)
{
  point_cloud_t* domain;
  stencil_t* neighborhoods;
  real_t* smoothing_lengths;
  static real_t h_over_dx[] = {1.0, 1.2, 1.5, 2.0, 2.4};
  make_lattice(10, 10, 10, h_over_dx[p], &domain, &neighborhoods, &smoothing_lengths);

  shape_function_kernel_t* W = simple_shape_function_kernel_new();
  shape_function_t* phi = mls_shape_function_new(p, W, domain, neighborhoods, smoothing_lengths);

  // Set up a constant field.
  int num_ghosts = stencil_num_ghosts(neighborhoods);
  int N = domain->num_points + num_ghosts;
  real_t one[N];
  for (int i = 0; i < N; ++i)
    one[i] = 1.0;

  // See if we can exactly reproduce it at random points within neighborhoods.
  rng_t* rng = host_rng_new();
  real_t dx = 0.1;
  for (int i = 0; i < domain->num_points; ++i)
  {
    shape_function_set_neighborhood(phi, i);
    point_t x = domain->points[i];
    bbox_t jitterbox = {.x1 = x.x - 0.5*dx, .x2 = x.x + 0.5*dx, 
                        .y1 = x.y - 0.5*dx, .y2 = x.y + 0.5*dx, 
                        .z1 = x.z - 0.5*dx, .z2 = x.z + 0.5*dx};
//    do
//    {
    point_randomize(&x, rng, &jitterbox);
//    }
//    while ((x.x < 0.5*dx) || (x.x > 1.0-0.5*dx) ||
//           (x.y < 0.5*dx) || (x.y > 1.0-0.5*dx) ||
//           (x.z < 0.5*dx) || (x.z > 1.0-0.5*dx));
           
    int N = shape_function_num_points(phi);
    real_t phi_val[N];
    vector_t phi_grad[N];
    shape_function_compute(phi, &x, phi_val, phi_grad);

    real_t val = 0.0; 
    vector_t grad = {.x = 0.0, .y = 0.0, .z = 0.0}; 
    int pos = 0, j, k = 0;
    while (stencil_next(neighborhoods, i, &pos, &j, NULL))
    {
      val += one[j] * phi_val[k];
      grad.x += one[j] * phi_grad[k].x;
      grad.y += one[j] * phi_grad[k].y;
      grad.z += one[j] * phi_grad[k].z;
      ++k;
    }

//printf("p = %d: %d points\n", p, N);
//printf("(%g, %g, %g), (%g, %g, %g), D/h = %g: %g %g %g %g\n", domain->points[i].x, domain->points[i].y, domain->points[i].z, x.x, x.y, x.z, point_distance(&domain->points[i], &x)/smoothing_lengths[i], val, grad.x, grad.y, grad.z);
    assert_true(fabs(val - 1.0) < 1e-14);
    assert_true(fabs(grad.x) < 1e-14);
    assert_true(fabs(grad.y) < 1e-14);
    assert_true(fabs(grad.z) < 1e-14);
  }

  // Clean up.
  shape_function_free(phi);
  point_cloud_free(domain);
  stencil_free(neighborhoods);
  polymec_free(smoothing_lengths);
}

void test_mls_shape_function_zero_consistency(void** state)
{
  for (int p = 0; p < 4; ++p)
    test_mls_shape_function_zero_consistency_p(state, p);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_mls_shape_function_ctors),
    unit_test(test_mls_shape_function_zero_consistency)
  };
  return run_tests(tests);
}
