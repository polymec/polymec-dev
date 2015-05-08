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
#include "core/unordered_set.h"
#include "geometry/cubic_lattice.h"
#include "geometry/create_uniform_mesh.h"
#include "model/mesh_stencils.h"

void test_NXxNYxNZ_stencil(void** state, 
                           MPI_Comm comm,
                           stencil_t* (*stencil_ctor)(mesh_t*, int),
                           int radius,
                           int nx, int ny, int nz,
                           int num_interior_neighbors, 
                           int num_boundary_neighbors, 
                           int num_edge_neighbors,
                           int num_corner_neighbors)
{
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  mesh_t* mesh = create_uniform_mesh(comm, nx, ny, nz, &bbox);
  cubic_lattice_t* lattice = cubic_lattice_new(nx, ny, nz);
  stencil_t* stencil = stencil_ctor(mesh, radius);
  for (int i = 0; i < nx; ++i)
  {
    for (int j = 0; j < ny; ++j)
    {
      for (int k = 0; k < nz; ++k)
      {
        // The cell score is a way of telling whether we're in the 
        // interior (0), on a boundary but not an edge/corner (1), on an 
        // edge but not a corner (2), or on a corner (3).
        int cell_score = 0;
        if ((i == 0) || (i == nx-1))
          ++cell_score;
        if ((j == 0) || (j == ny-1))
          ++cell_score;
        if ((k == 0) || (k == nz-1))
          ++cell_score;

        // Count the number of neighbors.
        int cell = cubic_lattice_cell(lattice, i, j, k);
        int num_neighbors = stencil_size(stencil, cell);
        if (cell_score == 0)
          assert_int_equal(num_interior_neighbors, num_neighbors);
        else if (cell_score == 1)
          assert_int_equal(num_boundary_neighbors, num_neighbors);
        else if (cell_score == 2)
          assert_int_equal(num_edge_neighbors, num_neighbors);
        else if (cell_score == 3)
          assert_int_equal(num_corner_neighbors, num_neighbors);

        // Make sure the neighbors are all valid, distinct cell indices.
        int_unordered_set_t* neighbors = int_unordered_set_new();
        int pos = 0, n;
        real_t weight;
        while (stencil_next(stencil, cell, &pos, &n, &weight))
        {
          assert_true(n >= 0);
          assert_true(n < mesh->num_cells);
          assert_false(int_unordered_set_contains(neighbors, n));
          int_unordered_set_insert(neighbors, n);
        }
        int_unordered_set_free(neighbors);
      }
    }
  }
  stencil_free(stencil);
  mesh_free(mesh);
}

void test_serial_1x1x1_cell_star_stencil(void** state)
{
  test_NXxNYxNZ_stencil(state, MPI_COMM_SELF, cell_star_stencil_new, 1, 1, 1, 1, 0, 0, 0, 0);
}

void test_serial_10x1x1_cell_star_stencil(void** state)
{
  test_NXxNYxNZ_stencil(state, MPI_COMM_SELF, cell_star_stencil_new, 1, 10, 1, 1, 0, 0, 2, 1);
}

void test_serial_10x10x1_cell_star_stencil(void** state)
{
  test_NXxNYxNZ_stencil(state, MPI_COMM_SELF, cell_star_stencil_new, 1, 10, 10, 1, 4, 4, 3, 2);
}

void test_serial_10x10x10_cell_star_stencil(void** state)
{
  test_NXxNYxNZ_stencil(state, MPI_COMM_SELF, cell_star_stencil_new, 1, 10, 10, 10, 6, 5, 4, 3);
}

void test_serial_1x1x1_cell_halo_stencil(void** state)
{
  test_NXxNYxNZ_stencil(state, MPI_COMM_SELF, cell_halo_stencil_new, 1, 1, 1, 1, 0, 0, 0, 0);
}

void test_serial_10x1x1_cell_halo_stencil(void** state)
{
  test_NXxNYxNZ_stencil(state, MPI_COMM_SELF, cell_halo_stencil_new, 1, 10, 1, 1, 0, 0, 2, 1);
}

void test_serial_10x10x1_cell_halo_stencil(void** state)
{
  test_NXxNYxNZ_stencil(state, MPI_COMM_SELF, cell_halo_stencil_new, 1, 10, 10, 1, 8, 8, 5, 3);
}

void test_serial_10x10x10_cell_halo_stencil(void** state)
{
  test_NXxNYxNZ_stencil(state, MPI_COMM_SELF, cell_halo_stencil_new, 1, 10, 10, 10, 26, 17, 11, 7);
}

void test_parallel_1x1x1_cell_star_stencil(void** state)
{
  test_NXxNYxNZ_stencil(state, MPI_COMM_WORLD, cell_star_stencil_new, 1, 1, 1, 1, 0, 0, 0, 0);
}

void test_parallel_10x1x1_cell_star_stencil(void** state)
{
  test_NXxNYxNZ_stencil(state, MPI_COMM_WORLD, cell_star_stencil_new, 1, 10, 1, 1, 0, 0, 2, 1);
}

void test_parallel_10x10x1_cell_star_stencil(void** state)
{
  test_NXxNYxNZ_stencil(state, MPI_COMM_WORLD, cell_star_stencil_new, 1, 10, 10, 1, 4, 4, 3, 2);
}

void test_parallel_10x10x10_cell_star_stencil(void** state)
{
  test_NXxNYxNZ_stencil(state, MPI_COMM_WORLD, cell_star_stencil_new, 1, 10, 10, 10, 6, 5, 4, 3);
}

void test_parallel_1x1x1_cell_halo_stencil(void** state)
{
  test_NXxNYxNZ_stencil(state, MPI_COMM_WORLD, cell_halo_stencil_new, 1, 1, 1, 1, 0, 0, 0, 0);
}

void test_parallel_10x1x1_cell_halo_stencil(void** state)
{
  test_NXxNYxNZ_stencil(state, MPI_COMM_WORLD, cell_halo_stencil_new, 1, 10, 1, 1, 0, 0, 2, 1);
}

void test_parallel_10x10x1_cell_halo_stencil(void** state)
{
  test_NXxNYxNZ_stencil(state, MPI_COMM_WORLD, cell_halo_stencil_new, 1, 10, 10, 1, 8, 8, 5, 3);
}

void test_parallel_10x10x10_cell_halo_stencil(void** state)
{
  test_NXxNYxNZ_stencil(state, MPI_COMM_WORLD, cell_halo_stencil_new, 1, 10, 10, 10, 26, 17, 11, 7);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_serial_1x1x1_cell_star_stencil),
    unit_test(test_serial_10x1x1_cell_star_stencil),
    unit_test(test_serial_10x10x1_cell_star_stencil),
    unit_test(test_serial_10x10x10_cell_star_stencil),
    unit_test(test_serial_1x1x1_cell_halo_stencil),
    unit_test(test_serial_10x1x1_cell_halo_stencil),
    unit_test(test_serial_10x10x1_cell_halo_stencil),
    unit_test(test_serial_10x10x10_cell_halo_stencil)
#if POLYMEC_HAVE_MPI
   ,unit_test(test_parallel_1x1x1_cell_star_stencil),
    unit_test(test_parallel_10x1x1_cell_star_stencil),
    unit_test(test_parallel_10x10x1_cell_star_stencil),
    unit_test(test_parallel_10x10x10_cell_star_stencil),
    unit_test(test_parallel_1x1x1_cell_halo_stencil),
    unit_test(test_parallel_10x1x1_cell_halo_stencil),
    unit_test(test_parallel_10x10x1_cell_halo_stencil),
    unit_test(test_parallel_10x10x10_cell_halo_stencil)
#endif
  };
  return run_tests(tests);
}
