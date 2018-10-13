// Copyright (c) 2012-2018, Jeffrey N. Johnson
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
#include "geometry/create_hex_prismesh.h"

static void test_create_hex_prismesh(void** state)
{
  // Create a uniform hex prismesh.
  MPI_Comm comm = MPI_COMM_WORLD;
  real_t h = 0.1, z1 = 0.0, z2 = 1.0;
  size_t radius = 5, nz = 10;
  prismesh_t* mesh = create_hex_prismesh(comm, radius, h,
                                         nz, z1, z2, false);
  real_t z1_, z2_;
  bool periodic;
  prismesh_get_z_info(mesh, &z1_, &z2_, &periodic);
  assert_true(reals_equal(z1_, z1));
  assert_true(reals_equal(z2_, z2));
  assert_false(periodic);
  assert_true(prismesh_comm(mesh) == comm);

  // Verify the mesh's topology.
  assert_true(prismesh_verify_topology(mesh, polymec_error));

  // Now check its connectivity.
  int local_num_cells = 0, local_num_faces = 0, 
      local_num_edges = 0, local_num_nodes = 0;
  prismesh_chunk_t* chunk;
  int pos = 0, xy_index, z_index;
  while (prismesh_next_chunk(mesh, &pos, &xy_index, &z_index, &chunk))
  {
    local_num_cells += (int)(chunk->num_columns * chunk->num_z_cells);
    local_num_faces += (int)(chunk->num_xy_faces * chunk->num_z_cells + 
                             (chunk->num_z_cells + 1));
    local_num_edges += (int)(chunk->num_xy_edges * (2 * chunk->num_z_cells + 1));
    local_num_nodes += (int)(chunk->num_xy_nodes * (chunk->num_z_cells + 1));
  }

  int global_num_cells, global_num_faces, global_num_edges, global_num_nodes;
  MPI_Allreduce(&local_num_cells, &global_num_cells, 1, MPI_INT, MPI_SUM, comm);
  MPI_Allreduce(&local_num_faces, &global_num_faces, 1, MPI_INT, MPI_SUM, comm);
  MPI_Allreduce(&local_num_edges, &global_num_edges, 1, MPI_INT, MPI_SUM, comm);
  MPI_Allreduce(&local_num_nodes, &global_num_nodes, 1, MPI_INT, MPI_SUM, comm);
  int nproc;
  MPI_Comm_size(comm, &nproc);
  assert_int_equal(10*10*10, global_num_cells);
  if (nproc > 1)
  {
    assert_true(global_num_faces > 3*10*10*11);
    assert_true(global_num_edges > 3*10*11*11);
    assert_true(global_num_nodes > 11*11*11);
  }
  else
  {
    assert_int_equal(3*10*10*11, global_num_faces);
    assert_int_equal(3*10*11*11, global_num_edges);
    assert_int_equal(11*11*11, global_num_nodes);
  }

  prismesh_free(mesh);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_create_hex_prismesh)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
