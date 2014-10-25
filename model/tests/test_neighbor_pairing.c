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
#include "core/kd_tree.h"
#include "core/array_utils.h"
#include "geometry/create_point_lattice.h"
#include "model/point_weight_function.h"
#include "model/neighbor_pairing.h"

// These functions implement a hat-shaped weight function with constant 
// spatial extent h.
static real_t hat_W(void* context, vector_t* y)
{
  real_t h = *((real_t*)context); // Spatial extent of W.
  real_t D = vector_mag(y);
  return (D < h) ? 1.0 - D/h : 0.0;
}

static vector_t hat_grad(void* context, vector_t* y)
{
  real_t h = *((real_t*)context); // Spatial extent of W.
  real_t D = vector_mag(y);
  static vector_t zero = {.x = 0.0, .y = 0.0, .z = 0.0};
  vector_t g = {.x = -y->x/(D*h), .y = -y->y/(D*h), .z = -y->z/(D*h)};
  return (D < h) ? g : zero;
}

static void hat_dtor(void* context)
{
  polymec_free(context);
}

static point_weight_function_t* hat_function_new(real_t h)
{
  real_t* Wh = polymec_malloc(sizeof(real_t));
  *Wh = h;
  point_weight_function_vtable vtable = {.value = hat_W, 
                                         .gradient = hat_grad, 
                                         .dtor = hat_dtor};
  return point_weight_function_new("Hat function", Wh, vtable);
}

static neighbor_pairing_t* create_simple_pairing(point_cloud_t* cloud,
                                                 real_t h)
{
  point_weight_function_t* W = hat_function_new(h);

  // Toss all the points into a kd-tree.
  kd_tree_t* tree = kd_tree_new(cloud->points, cloud->num_points);

  // Now do a neighbor search -- everything that falls within h of a point 
  // is a neighbor.
  int_array_t* pairs = int_array_new();
  real_array_t* weights = real_array_new();
  for (int i = 0; i < cloud->num_points; ++i)
  {
    point_t* xi = &cloud->points[i];
    int_slist_t* neighbors = kd_tree_within_radius(tree, &cloud->points[i], h);
    int_slist_node_t* n = NULL;
    int j;
    while (int_slist_next(neighbors, &n, &j))
    {
      if (j > i)
      {
        int_array_append(pairs, i);
        int_array_append(pairs, j);
        point_t* xj = &cloud->points[j];
        vector_t y;
        point_displacement(xj, xi, &y);
        real_array_append(weights, point_weight_function_value(W, &y));
      }
    }
    int_slist_free(neighbors);
  }
  kd_tree_free(tree);
  point_weight_function_free(W);

  exchanger_t* ex = exchanger_new(MPI_COMM_SELF);
  neighbor_pairing_t* pairing = neighbor_pairing_new("simple pairing", 
                                                     pairs->size/2,
                                                     pairs->data,
                                                     weights->data,
                                                     ex);
  int_array_release_data_and_free(pairs);
  real_array_release_data_and_free(weights);
  return pairing;
}

void test_point_lattice(void** state, 
                        MPI_Comm comm,
                        int nx, int ny, int nz, real_t h,
                        int num_interior_neighbors, 
                        int num_boundary_neighbors, 
                        int num_edge_neighbors,
                        int num_corner_neighbors)
{
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  point_cloud_t* cloud = create_uniform_point_lattice(comm, nx, ny, nz, &bbox);
  neighbor_pairing_t* pairing = create_simple_pairing(cloud, h);

  // Find the numbers of neighbors of each of the points. Do it collecting 
  // weights first, and then ignoring weights, and make sure we get the same
  // result.
  int num_neighbors1[nx*ny*nz], num_neighbors2[nx*ny*nz]; 
  memset(num_neighbors1, 0, nx*ny*nz * sizeof(int));
  memset(num_neighbors2, 0, nx*ny*nz * sizeof(int));
  int pos = 0, i, j;
  real_t wij;
  while (neighbor_pairing_next(pairing, &pos, &i, &j, &wij))
  {
    ++num_neighbors1[i];
    ++num_neighbors1[j];
  }
  pos = 0;
  while (neighbor_pairing_next(pairing, &pos, &i, &j, NULL))
  {
    ++num_neighbors2[i];
    ++num_neighbors2[j];
  }
  for (int i = 0; i < nx*ny*nz; ++i)
  {
    assert_int_equal(num_neighbors1[i], num_neighbors2[i]);
  }

  // Make bins of numbers of neighbors. There shouldn't be more than 4 
  // bins, and they should be (in ascending order): num_corner_neighbors, 
  // num_edge_neighbors, num_boundary_neighbors, num_interior_neighbors.
  int bins[1000]; // Up to 1000 neighbors (ridiculous).
  memset(bins, 0, 1000 * sizeof(int));
  for (int i = 0; i < nx*ny*nz; ++i)
    ++bins[num_neighbors1[i]];
  int num_nonempty_bins = 0;
  for (int i = 0; i < 1000; ++i)
  {
    if (bins[i] > 0)
      ++num_nonempty_bins;
  }
  assert_true(num_nonempty_bins <= 4);

  // Now order the bin counts.
  int bin_counts[4] = {-1, -1, -1, -1}, counter = 0;
  for (int i = 0; i < 1000; ++i)
  {
    if (bins[i] > 0)
      bin_counts[counter++] = i;
  }
  int_qsort(bin_counts, num_nonempty_bins);

  // Now check the neighbor counts against the reference ones we're given.
  assert_int_equal(num_corner_neighbors, bin_counts[0]);
  assert_true((num_edge_neighbors == bin_counts[0]) || 
              (num_edge_neighbors == bin_counts[1]));
  assert_true((num_boundary_neighbors == bin_counts[0]) || 
              (num_boundary_neighbors == bin_counts[1]) || 
              (num_boundary_neighbors == bin_counts[2]));
  assert_true((num_interior_neighbors == bin_counts[0]) || 
              (num_interior_neighbors == bin_counts[1]) || 
              (num_interior_neighbors == bin_counts[2]) ||
              (num_interior_neighbors == bin_counts[3]));

  // Clean up.
  neighbor_pairing_free(pairing);
  point_cloud_free(cloud);
}

void test_serial_1x1x1_lattice(void** state)
{
  test_point_lattice(state, MPI_COMM_SELF, 1, 1, 1, 0.1, 0, 0, 0, 0);
}

void test_serial_10x1x1_lattice(void** state)
{
  test_point_lattice(state, MPI_COMM_SELF, 10, 1, 1, 0.15, 2, 2, 2, 1);
}

void test_serial_10x10x1_lattice(void** state)
{
  test_point_lattice(state, MPI_COMM_SELF, 10, 10, 1, 0.15, 8, 8, 5, 3);
}

void test_serial_10x10x10_lattice(void** state)
{
  test_point_lattice(state, MPI_COMM_SELF, 10, 10, 10, 0.15, 18, 13, 9, 6);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_serial_1x1x1_lattice),
    unit_test(test_serial_10x1x1_lattice),
    unit_test(test_serial_10x10x1_lattice),
    unit_test(test_serial_10x10x10_lattice)
  };
  return run_tests(tests);
}
