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
#include "geometry/create_uniform_polymesh.h"

static void test_create_uniform_mesh(void** state)
{
  // Create a 10x10x10 uniform mesh.
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  polymesh_t* mesh = create_uniform_polymesh(MPI_COMM_WORLD, 10, 10, 10, &bbox);

  // Verify the mesh's topology.
  assert_true(polymesh_is_valid(mesh, NULL));

  // Now check its connectivity.
  int nproc;
  MPI_Comm_size(mesh->comm, &nproc);
  int num_cells, num_faces, num_edges, num_nodes;
  MPI_Allreduce(&mesh->num_cells, &num_cells, 1, MPI_INT, MPI_SUM, mesh->comm);
  MPI_Allreduce(&mesh->num_faces, &num_faces, 1, MPI_INT, MPI_SUM, mesh->comm);
  MPI_Allreduce(&mesh->num_edges, &num_edges, 1, MPI_INT, MPI_SUM, mesh->comm);
  MPI_Allreduce(&mesh->num_nodes, &num_nodes, 1, MPI_INT, MPI_SUM, mesh->comm);
  assert_int_equal(10*10*10, num_cells);
  if (nproc > 1)
  {
    assert_true(mesh->num_ghost_cells > 0);
    assert_true(num_faces > 3*10*10*11);
    assert_true(num_edges > 3*10*11*11);
    assert_true(num_nodes > 11*11*11);
  }
  else
  {
    assert_int_equal(0, mesh->num_ghost_cells);
    assert_int_equal(3*10*10*11, num_faces);
    assert_int_equal(3*10*11*11, num_edges);
    assert_int_equal(11*11*11, num_nodes);
  }

  polymesh_free(mesh);
}

static void test_create_uniform_mesh_on_rank(void** state)
{
  // Create a 10x10x10 rectilinear mesh on rank 0.
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  polymesh_t* mesh = create_uniform_polymesh_on_rank(MPI_COMM_WORLD, 0, 10, 10, 10, &bbox);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0)
  {
    assert_true(polymesh_is_valid(mesh, NULL));
    polymesh_free(mesh);
  }
  else
  {
    assert_true(mesh == NULL);
  }
}

static void test_plot_uniform_mesh_with_num_files(void** state, int num_files)
{
  // Create a uniform mesh.
  int nx = 10, ny = 10, nz = 10;
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  polymesh_t* mesh = create_uniform_polymesh(MPI_COMM_WORLD, nx, ny, nz, &bbox);

  // Clean up.
  polymesh_free(mesh);
}

static void test_plot_uniform_mesh_to_single_file(void** state)
{
  test_plot_uniform_mesh_with_num_files(state, 1);
}

static void test_plot_uniform_mesh_to_n_files(void** state)
{
  int nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  test_plot_uniform_mesh_with_num_files(state, nprocs);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_create_uniform_mesh),
    cmocka_unit_test(test_create_uniform_mesh_on_rank),
    cmocka_unit_test(test_plot_uniform_mesh_to_single_file),
    cmocka_unit_test(test_plot_uniform_mesh_to_n_files)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
