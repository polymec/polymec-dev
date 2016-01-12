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
#include "cmockery.h"
#include "core/array_utils.h"
#include "geometry/create_point_lattice.h"
#include "model/point_weight_function.h"
#include "model/neighbor_pairing.h"

// This creates a neighbor pairing using a hat function.
extern neighbor_pairing_t* create_simple_pairing(point_cloud_t* cloud, real_t h);

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
