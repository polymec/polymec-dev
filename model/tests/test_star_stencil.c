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
#include "core/unordered_set.h"
#include "core/timer.h"
#include "geometry/cubic_lattice.h"
#include "geometry/create_uniform_polymesh.h"
#include "model/polymesh_stencils.h"

static void check_stencil(void** state, 
                          polymesh_t* mesh,
                          int nx, int ny, int nz,
                          int num_interior_neighbors, 
                          int num_boundary_neighbors,
                          int num_edge_neighbors, 
                          int num_corner_neighbors,
                          stencil_t* stencil)
{
  START_FUNCTION_TIMER();

  MPI_Comm comm = mesh->comm;
  int rank, nprocs;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &nprocs);
  int num_cells[nprocs];
  MPI_Allgather(&mesh->num_cells, 1, MPI_INT, num_cells, 1, MPI_INT, comm);
  int cell_offset = 0;
  for (int p = 0; p < rank; ++p)
    cell_offset += num_cells[p];
  cubic_lattice_t* lattice = cubic_lattice_new(nx, ny, nz);

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
    size_t num_neighbors = stencil_size(stencil, c);
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
    while (stencil_next(stencil, c, &pos, &n))
    {
      assert_true(n >= 0);
      assert_false(int_unordered_set_contains(neighbors, n));
      int_unordered_set_insert(neighbors, n);
    }
    int_unordered_set_free(neighbors);
  }
  STOP_FUNCTION_TIMER();
}

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
  polymesh_t* mesh = create_uniform_polymesh(comm, nx, ny, nz, &bbox);
  stencil_t* stencil = cell_star_stencil_new(mesh, radius);
  check_stencil(state, mesh, nx, ny, nz, 
                num_interior_neighbors, num_boundary_neighbors,
                num_edge_neighbors, num_corner_neighbors,
                stencil);

  // Make a clone of this stencil and test it.
  stencil_t* stencil1 = stencil_clone(stencil);
  check_stencil(state, mesh, nx, ny, nz, 
                num_interior_neighbors, num_boundary_neighbors,
                num_edge_neighbors, num_corner_neighbors,
                stencil1);
  stencil_free(stencil1);

  // Serialize and deserialize the stencil.
  serializer_t* S = stencil_serializer();
  byte_array_t* bytes = byte_array_new();
  size_t offset = 0;
  serializer_write(S, stencil, bytes, &offset);
  offset = 0;
  stencil1 = serializer_read(S, bytes, &offset);
  check_stencil(state, mesh, nx, ny, nz, 
                num_interior_neighbors, num_boundary_neighbors,
                num_edge_neighbors, num_corner_neighbors,
                stencil1);
  stencil_free(stencil1);
  byte_array_free(bytes);

  // Create a graph from the stencil.
  adj_graph_t* G = stencil_as_graph(stencil);
  assert_int_equal(adj_graph_num_vertices(G), stencil_num_indices(stencil));
  adj_graph_free(G);

  // Create a matrix sparsity from the stencil.
  matrix_sparsity_t* sp = matrix_sparsity_from_stencil(stencil);
  assert_int_equal(matrix_sparsity_num_local_rows(sp), stencil_num_indices(stencil));
  matrix_sparsity_free(sp);

  // Clean up.
  stencil_free(stencil);
  polymesh_free(mesh);
}

static void test_serial_1x1x1_cell_star_stencil(void** state)
{
  START_FUNCTION_TIMER();
  test_NXxNYxNZ_star_stencil(state, MPI_COMM_SELF, 1, 1, 1, 1, 0, 0, 0, 0);
  STOP_FUNCTION_TIMER();
}

static void test_serial_10x1x1_cell_star_stencil(void** state)
{
  START_FUNCTION_TIMER();
  test_NXxNYxNZ_star_stencil(state, MPI_COMM_SELF, 1, 10, 1, 1, 0, 0, 2, 1);
  STOP_FUNCTION_TIMER();
}

static void test_serial_10x10x1_cell_star_stencil(void** state)
{
  START_FUNCTION_TIMER();
  test_NXxNYxNZ_star_stencil(state, MPI_COMM_SELF, 1, 10, 10, 1, 4, 4, 3, 2);
  STOP_FUNCTION_TIMER();
}

static void test_serial_10x10x10_cell_star_stencil(void** state)
{
  test_NXxNYxNZ_star_stencil(state, MPI_COMM_SELF, 1, 10, 10, 10, 6, 5, 4, 3);
  START_FUNCTION_TIMER();
  STOP_FUNCTION_TIMER();
}

#if POLYMEC_HAVE_MPI
static void test_parallel_10x1x1_cell_star_stencil(void** state)
{
  START_FUNCTION_TIMER();
  test_NXxNYxNZ_star_stencil(state, MPI_COMM_WORLD, 1, 10, 1, 1, 0, 0, 2, 1);
  STOP_FUNCTION_TIMER();
}

static void test_parallel_10x10x1_cell_star_stencil(void** state)
{
  START_FUNCTION_TIMER();
  test_NXxNYxNZ_star_stencil(state, MPI_COMM_WORLD, 1, 10, 10, 1, 4, 4, 3, 2);
  STOP_FUNCTION_TIMER();
}

static void test_parallel_10x10x10_cell_star_stencil(void** state)
{
  START_FUNCTION_TIMER();
  test_NXxNYxNZ_star_stencil(state, MPI_COMM_WORLD, 1, 10, 10, 10, 6, 5, 4, 3);
  STOP_FUNCTION_TIMER();
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
