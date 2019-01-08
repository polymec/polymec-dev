// Copyright (c) 2012-2019, Jeffrey N. Johnson
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
#include "core/array_utils.h"
#include "geometry/partition_point_cloud.h"
#include "geometry/create_point_lattice.h"
#include "model/neighbor_pairing.h"

// This creates a neighbor pairing using a hat function.
extern neighbor_pairing_t* create_simple_pairing(point_cloud_t* cloud, real_t h);

static void test_serial_point_lattice(void** state, 
                                      int nx, int ny, int nz, real_t h,
                                      int num_interior_neighbors, 
                                      int num_boundary_neighbors, 
                                      int num_edge_neighbors,
                                      int num_corner_neighbors)
{
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  point_cloud_t* cloud = create_uniform_point_lattice(MPI_COMM_SELF, nx, ny, nz, &bbox);
  neighbor_pairing_t* pairing = create_simple_pairing(cloud, h);

  // Find the numbers of neighbors of each of the points.
  int num_neighbors[nx*ny*nz];
  memset(num_neighbors, 0, nx*ny*nz * sizeof(int));
  {
    int pos = 0, i, j;
    while (neighbor_pairing_next(pairing, &pos, &i, &j))
    {
      ++num_neighbors[i];
      ++num_neighbors[j];
    }
  }
  for (int i = 0; i < nx*ny*nz; ++i)
  {
    assert_true(num_neighbors[i] >= 0);
  }

  // Make bins of numbers of neighbors. There shouldn't be more than 4 
  // bins, and they should be (in ascending order): num_corner_neighbors, 
  // num_edge_neighbors, num_boundary_neighbors, num_interior_neighbors.
  int bins[1000]; // Up to 1000 neighbors (ridiculous).
  memset(bins, 0, 1000 * sizeof(int));
  for (int i = 0; i < nx*ny*nz; ++i)
    ++bins[num_neighbors[i]];
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

static void test_parallel_point_lattice(void** state,
                                        int nx, int ny, int nz, real_t R0,
                                        int num_interior_neighbors, 
                                        int num_boundary_neighbors, 
                                        int num_edge_neighbors,
                                        int num_corner_neighbors)
{
  // Set up and partition the point cloud.
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  point_cloud_t* cloud = create_uniform_point_lattice(MPI_COMM_SELF, nx, ny, nz, &bbox);
  assert_true(partition_point_cloud(&cloud, MPI_COMM_WORLD, NULL, 0.10, NULL, 0));

  // Now set up the point radius field and use it to find neighbor pairs.
  int num_ghosts;
  real_t R[cloud->num_points];
  for (int i = 0; i < cloud->num_points + cloud->num_ghosts; ++i)
    R[i] = R0;
  neighbor_pairing_t* pairing = distance_based_neighbor_pairing_new(cloud, R, &num_ghosts);
  point_cloud_set_num_ghosts(cloud, num_ghosts);

  // Find the numbers of neighbors of each of the points. 
  size_t num_points = cloud->num_points;
  int num_neighbors[num_points];
  memset(num_neighbors, 0, num_points * sizeof(int));
  {
    int pos = 0, i, j;
    while (neighbor_pairing_next(pairing, &pos, &i, &j))
    {
      ++num_neighbors[i];
      if (j < num_points)
        ++num_neighbors[j];
    }
  }
  for (int i = 0; i < num_points; ++i)
  {
    assert_true(num_neighbors[i] >= 0); 
  }

  // Make bins of numbers of neighbors. There shouldn't be more than 4 
  // bins, and they should be (in ascending order): num_corner_neighbors, 
  // num_edge_neighbors, num_boundary_neighbors, num_interior_neighbors.
  int bins[1000]; // Up to 1000 neighbors (ridiculous).
  memset(bins, 0, 1000 * sizeof(int));
  for (int i = 0; i < num_points; ++i)
    ++bins[num_neighbors[i]];
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

  // Before we go, create a stencil from this pairing.
  stencil_t* S = stencil_from_point_cloud_and_neighbors(cloud, pairing);

  // Convert back to a pairing.
  neighbor_pairing_t* pairing1 = neighbor_pairing_from_stencil(S);
  assert_int_equal(strcmp(pairing1->name, pairing->name), 0);
  assert_int_equal(pairing1->num_pairs, pairing->num_pairs);
  for (int i = 0; i < pairing->num_pairs; ++i)
  {
    assert_int_equal(pairing1->pairs[2*i], pairing->pairs[2*i]);
    assert_int_equal(pairing1->pairs[2*i+1], pairing->pairs[2*i+1]);
  }

  // Clean up.
  neighbor_pairing_free(pairing1);
  stencil_free(S);
  neighbor_pairing_free(pairing);
  point_cloud_free(cloud);
}

static void test_serial_1x1x1_lattice(void** state)
{
  test_serial_point_lattice(state, 1, 1, 1, 0.1, 0, 0, 0, 0);
}

static void test_serial_10x1x1_lattice(void** state)
{
  test_serial_point_lattice(state, 10, 1, 1, 0.15, 2, 2, 2, 1);
}

static void test_serial_10x10x1_lattice(void** state)
{
  test_serial_point_lattice(state, 10, 10, 1, 0.15, 8, 8, 5, 3);
}

static void test_serial_10x10x10_lattice(void** state)
{
  test_serial_point_lattice(state, 10, 10, 10, 0.15, 18, 13, 9, 6);
}

static void test_parallel_10x10x1_lattice(void** state)
{
  test_parallel_point_lattice(state, 10, 10, 1, 0.15, 8, 8, 5, 3);
}

//static void test_parallel_10x10x10_lattice(void** state)
//{
//  test_parallel_point_lattice(state, 10, 10, 10, 0.15, 18, 13, 9, 6);
//}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_serial_1x1x1_lattice),
    cmocka_unit_test(test_serial_10x1x1_lattice),
    cmocka_unit_test(test_serial_10x10x1_lattice),
    cmocka_unit_test(test_serial_10x10x10_lattice),
    cmocka_unit_test(test_parallel_10x10x1_lattice),
    cmocka_unit_test(test_parallel_10x10x1_lattice),
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
