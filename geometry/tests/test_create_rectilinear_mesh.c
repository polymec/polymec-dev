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
#include "core/array_utils.h"
#include "core/silo_file.h"
#include "geometry/create_rectilinear_mesh.h"

void test_create_rectilinear_mesh(void** state)
{
  // Create a 10x10x10 rectilinear mesh.
  double xs[] = {0.0, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0, 256.0, 512.0};
  double ys[] = {0.0, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0, 256.0, 512.0};
  double zs[] = {0.0, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0, 256.0, 512.0};
  mesh_t* mesh = create_rectilinear_mesh(MPI_COMM_WORLD, xs, 11, ys, 11, zs, 11);
  assert_true(mesh_verify_topology(mesh, polymec_error));

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

  mesh_free(mesh);
}

void test_plot_rectilinear_mesh(void** state)
{
  // Create a 4x4x4 rectilinear mesh.
  double xs[] = {0.0, 1.0, 2.0, 4.0, 8.0};
  double ys[] = {0.0, 1.0, 2.0, 4.0, 8.0};
  double zs[] = {0.0, 1.0, 2.0, 4.0, 8.0};
  mesh_t* mesh = create_rectilinear_mesh(MPI_COMM_WORLD, xs, 5, ys, 5, zs, 5);

  // Plot it.
  double ones[4*4*4];
  for (int c = 0; c < 4*4*4; ++c)
    ones[c] = 1.0*c;
  silo_file_t* silo = silo_file_new(MPI_COMM_WORLD, "rectilinear_4x4x4", "", 1, 0, 0, 0.0);
  silo_file_write_mesh(silo, "mesh", mesh);
  silo_file_write_scalar_cell_field(silo, "solution", "mesh", ones, NULL);
  silo_file_close(silo);

  // Clean up.
  mesh_free(mesh);
}

static void check_cell_face_connectivity(void** state, 
                                         mesh_t* mesh,
                                         index_t* global_cell_indices, 
                                         index_t global_cell_index, 
                                         int cell_face_index, 
                                         index_t global_opp_cell_index)
{
  ASSERT(cell_face_index >= 0);
  ASSERT(cell_face_index < 6);

  // Find the given global cell index within our array.
  int cell_array_len = mesh->num_cells + mesh->num_ghost_cells;
  index_t* cell_p = index_lsearch(global_cell_indices, cell_array_len, global_cell_index);
  assert_true(cell_p != NULL);
  int cell = cell_p - global_cell_indices;
  int face = mesh->cell_faces[6*cell + cell_face_index];
  if (face < 0) 
    face = ~face;
  int opp_cell = mesh_face_opp_cell(mesh, face, cell);
  if ((opp_cell == -1) && (global_opp_cell_index == (index_t)(-1)))
    assert_int_equal(global_opp_cell_index, opp_cell);
  else
    assert_int_equal(global_opp_cell_index, global_cell_indices[opp_cell]); 
}

void test_3proc_4x4x1_mesh(void** state)
{
  // Create a 4x4x1 rectilinear mesh, which isn't a big deal in and of itself, 
  // but does seem to exhibit problems on certain numbers of processes.
  double xs[] = {0.0, 1.0, 2.0, 3.0, 4.0};
  double ys[] = {0.0, 1.0, 2.0, 3.0, 4.0};
  double zs[] = {0.0, 1.0};
  mesh_t* mesh = create_rectilinear_mesh(MPI_COMM_WORLD, xs, 5, ys, 5, zs, 2);
  assert_true(mesh_verify_topology(mesh, polymec_error));

  int rank;
  MPI_Comm_rank(mesh->comm, &rank);

  // Create global cell indices.
  index_t vtx_dist[4] = {0, 5, 10, 16};
  index_t G[mesh->num_cells + mesh->num_ghost_cells];
  for (int i = 0; i < mesh->num_cells; ++i)
    G[i] = vtx_dist[rank] + i;
  exchanger_exchange(mesh_exchanger(mesh), G, 1, 0, MPI_INDEX_T);

  if (rank == 0)
  {
    assert_int_equal(5, mesh->num_cells);
    assert_int_equal(5, mesh->num_ghost_cells);

    // Check the layouts of local and ghost cells.
    // Note that ghost cells are created by walking through the 
    // local cells and trying to add ghosts in each of the -/+ x/y/z 
    // directions (in that order).
    index_t G0[] = {0, 1, 2, 3, 4, 5, 5, 6, 7, 8};
    for (int i = 0; i < 10; ++i)
      assert_int_equal(G0[i], G[i]);

    // Check that cell->face and face->cell connectivity are 
    // correctly and consistently done.
    check_cell_face_connectivity(state, mesh, G, 0, 0, -1);
    check_cell_face_connectivity(state, mesh, G, 0, 1,  1);
    check_cell_face_connectivity(state, mesh, G, 0, 2, -1);
    check_cell_face_connectivity(state, mesh, G, 0, 3,  4);
    check_cell_face_connectivity(state, mesh, G, 0, 4, -1);
    check_cell_face_connectivity(state, mesh, G, 0, 5, -1);

    check_cell_face_connectivity(state, mesh, G, 1, 0,  0);
    check_cell_face_connectivity(state, mesh, G, 1, 1,  2);
    check_cell_face_connectivity(state, mesh, G, 1, 2, -1);
    check_cell_face_connectivity(state, mesh, G, 1, 3,  5);
    check_cell_face_connectivity(state, mesh, G, 1, 4, -1);
    check_cell_face_connectivity(state, mesh, G, 1, 5, -1);

    check_cell_face_connectivity(state, mesh, G, 2, 0,  1);
    check_cell_face_connectivity(state, mesh, G, 2, 1,  3);
    check_cell_face_connectivity(state, mesh, G, 2, 2, -1);
    check_cell_face_connectivity(state, mesh, G, 2, 3,  6);
    check_cell_face_connectivity(state, mesh, G, 2, 4, -1);
    check_cell_face_connectivity(state, mesh, G, 2, 5, -1);

    check_cell_face_connectivity(state, mesh, G, 3, 0,  2);
    check_cell_face_connectivity(state, mesh, G, 3, 1, -1);
    check_cell_face_connectivity(state, mesh, G, 3, 2, -1);
    check_cell_face_connectivity(state, mesh, G, 3, 3,  7);
    check_cell_face_connectivity(state, mesh, G, 3, 4, -1);
    check_cell_face_connectivity(state, mesh, G, 3, 5, -1);

    check_cell_face_connectivity(state, mesh, G, 4, 0, -1);
    check_cell_face_connectivity(state, mesh, G, 4, 1,  5);
    check_cell_face_connectivity(state, mesh, G, 4, 2,  0);
    check_cell_face_connectivity(state, mesh, G, 4, 3,  8);
    check_cell_face_connectivity(state, mesh, G, 4, 4, -1);
    check_cell_face_connectivity(state, mesh, G, 4, 5, -1);
  }
  else if (rank == 1)
  {
    assert_int_equal(5, mesh->num_cells);
    assert_int_equal(10, mesh->num_ghost_cells);

    // Check the layouts of local and ghost cells.
    // Note that ghost cells are created by walking through the 
    // local cells and trying to add ghosts in each of the -/+ x/y/z 
    // directions (in that order).
    index_t G1[] = {5, 6, 7, 8, 9, 1, 2, 3, 10, 4, 10, 4, 11, 12, 13};
    for (int i = 0; i < 15; ++i)
      assert_int_equal(G1[i], G[i]);

    // Check that cell->face and face->cell connectivity are 
    // correctly and consistently done.
    check_cell_face_connectivity(state, mesh, G, 5, 0,  4);
    check_cell_face_connectivity(state, mesh, G, 5, 1,  6);
    check_cell_face_connectivity(state, mesh, G, 5, 2,  1);
    check_cell_face_connectivity(state, mesh, G, 5, 3,  9);
    check_cell_face_connectivity(state, mesh, G, 5, 4, -1);
    check_cell_face_connectivity(state, mesh, G, 5, 5, -1);

    check_cell_face_connectivity(state, mesh, G, 6, 0,  5);
    check_cell_face_connectivity(state, mesh, G, 6, 1,  7);
    check_cell_face_connectivity(state, mesh, G, 6, 2,  2);
    check_cell_face_connectivity(state, mesh, G, 6, 3, 10);
    check_cell_face_connectivity(state, mesh, G, 6, 4, -1);
    check_cell_face_connectivity(state, mesh, G, 6, 5, -1);

    check_cell_face_connectivity(state, mesh, G, 7, 0,  6);
    check_cell_face_connectivity(state, mesh, G, 7, 1, -1);
    check_cell_face_connectivity(state, mesh, G, 7, 2,  3);
    check_cell_face_connectivity(state, mesh, G, 7, 3, 11);
    check_cell_face_connectivity(state, mesh, G, 7, 4, -1);
    check_cell_face_connectivity(state, mesh, G, 7, 5, -1);

    check_cell_face_connectivity(state, mesh, G, 8, 0, -1);
    check_cell_face_connectivity(state, mesh, G, 8, 1,  9);
    check_cell_face_connectivity(state, mesh, G, 8, 2,  4);
    check_cell_face_connectivity(state, mesh, G, 8, 3, 12);
    check_cell_face_connectivity(state, mesh, G, 8, 4, -1);
    check_cell_face_connectivity(state, mesh, G, 8, 5, -1);

    check_cell_face_connectivity(state, mesh, G, 9, 0,  8);
    check_cell_face_connectivity(state, mesh, G, 9, 1, 10);
    check_cell_face_connectivity(state, mesh, G, 9, 2,  5);
    check_cell_face_connectivity(state, mesh, G, 9, 3, 13);
    check_cell_face_connectivity(state, mesh, G, 9, 4, -1);
    check_cell_face_connectivity(state, mesh, G, 9, 5, -1);
  }
  else if (rank == 2)
  {
    assert_int_equal(6, mesh->num_cells);
    assert_int_equal(5, mesh->num_ghost_cells);

    // Check the layouts of local and ghost cells.
    // Note that ghost cells are created by walking through the 
    // local cells and trying to add ghosts in each of the -/+ x/y/z 
    // directions (in that order).
    index_t G2[] = {10, 11, 12, 13, 14, 15, 6, 7, 8, 9, 9};
    for (int i = 0; i < 11; ++i)
      assert_int_equal(G2[i], G[i]);

    // Check that cell->face and face->cell connectivity are 
    // correctly and consistently done.
    check_cell_face_connectivity(state, mesh, G, 10, 0,  9);
    check_cell_face_connectivity(state, mesh, G, 10, 1, 11);
    check_cell_face_connectivity(state, mesh, G, 10, 2,  6);
    check_cell_face_connectivity(state, mesh, G, 10, 3, 14);
    check_cell_face_connectivity(state, mesh, G, 10, 4, -1);
    check_cell_face_connectivity(state, mesh, G, 10, 5, -1);

    check_cell_face_connectivity(state, mesh, G, 11, 0, 10);
    check_cell_face_connectivity(state, mesh, G, 11, 1, -1);
    check_cell_face_connectivity(state, mesh, G, 11, 2,  7);
    check_cell_face_connectivity(state, mesh, G, 11, 3, 15);
    check_cell_face_connectivity(state, mesh, G, 11, 4, -1);
    check_cell_face_connectivity(state, mesh, G, 11, 5, -1);

    check_cell_face_connectivity(state, mesh, G, 12, 0, -1);
    check_cell_face_connectivity(state, mesh, G, 12, 1, 13);
    check_cell_face_connectivity(state, mesh, G, 12, 2,  8);
    check_cell_face_connectivity(state, mesh, G, 12, 3, -1);
    check_cell_face_connectivity(state, mesh, G, 12, 4, -1);
    check_cell_face_connectivity(state, mesh, G, 12, 5, -1);

    check_cell_face_connectivity(state, mesh, G, 13, 0, 12);
    check_cell_face_connectivity(state, mesh, G, 13, 1, 14);
    check_cell_face_connectivity(state, mesh, G, 13, 2,  9);
    check_cell_face_connectivity(state, mesh, G, 13, 3, -1);
    check_cell_face_connectivity(state, mesh, G, 13, 4, -1);
    check_cell_face_connectivity(state, mesh, G, 13, 5, -1);

    check_cell_face_connectivity(state, mesh, G, 14, 0, 13);
    check_cell_face_connectivity(state, mesh, G, 14, 1, 15);
    check_cell_face_connectivity(state, mesh, G, 14, 2, 10);
    check_cell_face_connectivity(state, mesh, G, 14, 3, -1);
    check_cell_face_connectivity(state, mesh, G, 14, 4, -1);
    check_cell_face_connectivity(state, mesh, G, 14, 5, -1);
  }

  mesh_free(mesh);
}

void test_problematic_meshes(void** state)
{
  int nproc;
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  if (nproc == 3)
  {
    // Check for correctness on the problematic 3-process 4x4x1 case.
    test_3proc_4x4x1_mesh(state);
  }
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_create_rectilinear_mesh),
    cmocka_unit_test(test_plot_rectilinear_mesh),
    cmocka_unit_test(test_problematic_meshes)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
