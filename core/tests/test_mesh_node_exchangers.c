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
#include "core/mesh.h"
#include "core/partition_mesh.h"
#include "geometry/create_uniform_mesh.h"

void test_nv_node_exchanger_on_line(void** state)
{
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  int nx = 8;
  mesh_t* mesh = create_uniform_mesh(MPI_COMM_SELF, nx, 1, 1, &bbox);
  partition_mesh(&mesh, MPI_COMM_WORLD, NULL, 0.0);
  int node_offsets[mesh->num_nodes+1];
  exchanger_t* ex = mesh_nv_node_exchanger_new(mesh, node_offsets);
exchanger_fprintf(ex, stdout);
  exchanger_verify(ex, polymec_error);

  int nprocs, rank;
  MPI_Comm_size(mesh->comm, &nprocs);
  MPI_Comm_rank(mesh->comm, &rank);

  if (nprocs == 2)
  {
    if (rank == 0)
    {
      assert_int_equal(1, exchanger_num_sends(ex));
      int* nodes, num_nodes;
      exchanger_get_send(ex, 1, &nodes, &num_nodes);
      assert_int_equal(4, num_nodes);
    }
    else
    {
      assert_int_equal(1, exchanger_num_receives(ex));
      int* nodes, num_nodes;
      exchanger_get_receive(ex, 0, &nodes, &num_nodes);
      assert_int_equal(4, num_nodes);
    }
  }

  exchanger_free(ex);
  mesh_free(mesh);
}

void test_1v_node_exchanger_on_line(void** state)
{
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  int nx = 8;
  mesh_t* mesh = create_uniform_mesh(MPI_COMM_SELF, nx, 1, 1, &bbox);
  partition_mesh(&mesh, MPI_COMM_WORLD, NULL, 0.0);
  exchanger_t* ex = mesh_1v_node_exchanger_new(mesh);

  int nprocs, rank;
  MPI_Comm_size(mesh->comm, &nprocs);
  MPI_Comm_rank(mesh->comm, &rank);

  if (nprocs == 2)
  {
    if (rank == 0)
    {
      assert_int_equal(1, exchanger_num_sends(ex));
      int* nodes, num_nodes;
      exchanger_get_send(ex, 1, &nodes, &num_nodes);
      assert_int_equal(4, num_nodes);
    }
    else
    {
      assert_int_equal(1, exchanger_num_receives(ex));
      int* nodes, num_nodes;
      exchanger_get_receive(ex, 0, &nodes, &num_nodes);
      assert_int_equal(4, num_nodes);
    }
  }

  exchanger_free(ex);
  mesh_free(mesh);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_nv_node_exchanger_on_line),
    unit_test(test_1v_node_exchanger_on_line)
  };
  return run_tests(tests);
}
