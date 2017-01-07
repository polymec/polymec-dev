// Copyright (c) 2012-2017, Jeffrey N. Johnson
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
#include "core/unordered_set.h"
#include "geometry/cubic_lattice.h"
#include "geometry/create_uniform_mesh.h"
#include "model/mesh_stencils.h"

static void test_NXxNYxNZ_star_stencil(void** state, 
                                       MPI_Comm comm,
                                       int radius,
                                       int nx, int ny, int nz,
                                       int num_interior_neighbors, 
                                       int num_boundary_neighbors, 
                                       int num_edge_neighbors,
                                       int num_corner_neighbors)
{
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  mesh_t* mesh = create_uniform_mesh(comm, nx, ny, nz, &bbox);
  int rank, nprocs;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &nprocs);
  int num_cells[nprocs];
  MPI_Allgather(&mesh->num_cells, 1, MPI_INT, num_cells, 1, MPI_INT, comm);
  int cell_offset = 0;
  for (int p = 0; p < rank; ++p)
    cell_offset += num_cells[p];
  cubic_lattice_t* lattice = cubic_lattice_new(nx, ny, nz);
  stencil_t* stencil = cell_star_stencil_new(mesh, radius);
  for (int cell = cell_offset; cell < cell_offset + mesh->num_cells; ++cell)
  {
    index_t i, j, k;
    cubic_lattice_get_cell_triple(lattice, (index_t)cell, &i, &j, &k);

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
    int c = cell - cell_offset;
    int num_neighbors = stencil_size(stencil, c);
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
    while (stencil_next(stencil, c, &pos, &n, &weight))
    {
      assert_true(n >= 0);
      assert_false(int_unordered_set_contains(neighbors, n));
      int_unordered_set_insert(neighbors, n);
    }
    int_unordered_set_free(neighbors);
  }
  stencil_free(stencil);
  mesh_free(mesh);
}

static void test_serial_1x1x1_cell_star_stencil(void** state)
{
  test_NXxNYxNZ_star_stencil(state, MPI_COMM_SELF, 1, 1, 1, 1, 0, 0, 0, 0);
}

static void test_serial_10x1x1_cell_star_stencil(void** state)
{
  test_NXxNYxNZ_star_stencil(state, MPI_COMM_SELF, 1, 10, 1, 1, 0, 0, 2, 1);
}

static void test_serial_10x10x1_cell_star_stencil(void** state)
{
  test_NXxNYxNZ_star_stencil(state, MPI_COMM_SELF, 1, 10, 10, 1, 4, 4, 3, 2);
}

static void test_serial_10x10x10_cell_star_stencil(void** state)
{
  test_NXxNYxNZ_star_stencil(state, MPI_COMM_SELF, 1, 10, 10, 10, 6, 5, 4, 3);
}

#if POLYMEC_HAVE_MPI
static void test_parallel_10x1x1_cell_star_stencil(void** state)
{
  test_NXxNYxNZ_star_stencil(state, MPI_COMM_WORLD, 1, 10, 1, 1, 0, 0, 2, 1);
}

static void test_parallel_10x10x1_cell_star_stencil(void** state)
{
  test_NXxNYxNZ_star_stencil(state, MPI_COMM_WORLD, 1, 10, 10, 1, 4, 4, 3, 2);
}

static void test_parallel_10x10x10_cell_star_stencil(void** state)
{
  test_NXxNYxNZ_star_stencil(state, MPI_COMM_WORLD, 1, 10, 10, 10, 6, 5, 4, 3);
}
#endif

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_serial_1x1x1_cell_star_stencil),
    cmocka_unit_test(test_serial_10x1x1_cell_star_stencil),
    cmocka_unit_test(test_serial_10x10x1_cell_star_stencil),
    cmocka_unit_test(test_serial_10x10x10_cell_star_stencil)
#if POLYMEC_HAVE_MPI
   ,cmocka_unit_test(test_parallel_10x1x1_cell_star_stencil),
    cmocka_unit_test(test_parallel_10x10x1_cell_star_stencil),
    cmocka_unit_test(test_parallel_10x10x10_cell_star_stencil)
#endif
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
