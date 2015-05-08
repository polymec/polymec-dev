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
  int nprocs, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  int nx = 8;
  mesh_t* mesh = create_uniform_mesh(MPI_COMM_SELF, nx, 1, 1, &bbox);
  partition_mesh(&mesh, MPI_COMM_WORLD, NULL, 0.0);
  int node_offsets[mesh->num_nodes+1];
  exchanger_t* ex = mesh_nv_node_exchanger_new(mesh, node_offsets);
  int node_owners[node_offsets[mesh->num_nodes]];
  for (int n = 0; n < node_offsets[mesh->num_nodes]; ++n)
    node_owners[n] = -1;
  for (int n = 0; n < mesh->num_nodes; ++n)
    node_owners[node_offsets[n]] = rank;
  exchanger_exchange(ex, node_owners, 1, 0, MPI_INT);
  for (int n = 0; n < mesh->num_nodes; ++n)
  {
    assert_int_equal(node_owners[node_offsets[n]], rank);
    for (int nn = node_offsets[n]+1; nn < node_offsets[n+1]; ++nn)
    {
      assert_true(node_owners[nn] >= 0);
      assert_true(node_owners[nn] < nprocs);
      assert_int_not_equal(node_owners[nn], rank);
    }
  }
  exchanger_fprintf(ex, stdout);
  exchanger_verify(ex, polymec_error);

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

  MPI_Barrier(MPI_COMM_WORLD);
}

void test_1v_node_exchanger_on_line(void** state)
{
  int nprocs, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  int nx = 8;
  mesh_t* mesh = create_uniform_mesh(MPI_COMM_SELF, nx, 1, 1, &bbox);
  partition_mesh(&mesh, MPI_COMM_WORLD, NULL, 0.0);
  exchanger_t* ex = mesh_1v_node_exchanger_new(mesh);
  exchanger_fprintf(ex, stdout);
  exchanger_verify(ex, polymec_error);

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

  MPI_Barrier(MPI_COMM_WORLD);
}

void test_nv_node_exchanger_in_plane(void** state)
{
  int nprocs, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  int nx = 8, ny = 8;
  mesh_t* mesh = create_uniform_mesh(MPI_COMM_SELF, nx, ny, 1, &bbox);
  partition_mesh(&mesh, MPI_COMM_WORLD, NULL, 0.0);
  int node_offsets[mesh->num_nodes+1];
  exchanger_t* ex = mesh_nv_node_exchanger_new(mesh, node_offsets);
  int node_owners[node_offsets[mesh->num_nodes]];
  for (int n = 0; n < node_offsets[mesh->num_nodes]; ++n)
    node_owners[n] = -1;
  for (int n = 0; n < mesh->num_nodes; ++n)
    node_owners[node_offsets[n]] = rank;
  exchanger_exchange(ex, node_owners, 1, 0, MPI_INT);
  for (int n = 0; n < mesh->num_nodes; ++n)
  {
    assert_int_equal(node_owners[node_offsets[n]], rank);
    for (int nn = node_offsets[n]+1; nn < node_offsets[n+1]; ++nn)
    {
      assert_true(node_owners[nn] >= 0);
      assert_true(node_owners[nn] < nprocs);
      assert_int_not_equal(node_owners[nn], rank);
    }
  }
  exchanger_fprintf(ex, stdout);
  exchanger_verify(ex, polymec_error);
  exchanger_free(ex);
  mesh_free(mesh);

  MPI_Barrier(MPI_COMM_WORLD);
}

void test_1v_node_exchanger_in_plane(void** state)
{
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  int nx = 4, ny = 4;
  mesh_t* mesh = create_uniform_mesh(MPI_COMM_SELF, nx, ny, 1, &bbox);
  partition_mesh(&mesh, MPI_COMM_WORLD, NULL, 0.0);
  exchanger_t* ex = mesh_1v_node_exchanger_new(mesh);
  exchanger_fprintf(ex, stdout);
  exchanger_verify(ex, polymec_error);
  exchanger_free(ex);
  mesh_free(mesh);

  MPI_Barrier(MPI_COMM_WORLD);
}

void test_nv_node_exchanger_in_cube(void** state)
{
  int nprocs, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  int nx = 4, ny = 4, nz = 4;
  mesh_t* mesh = create_uniform_mesh(MPI_COMM_SELF, nx, ny, nz, &bbox);
  partition_mesh(&mesh, MPI_COMM_WORLD, NULL, 0.0);
  int node_offsets[mesh->num_nodes+1];
  exchanger_t* ex = mesh_nv_node_exchanger_new(mesh, node_offsets);
  int node_owners[node_offsets[mesh->num_nodes]];
  for (int n = 0; n < node_offsets[mesh->num_nodes]; ++n)
    node_owners[n] = -1;
  for (int n = 0; n < mesh->num_nodes; ++n)
    node_owners[node_offsets[n]] = rank;
  exchanger_exchange(ex, node_owners, 1, 0, MPI_INT);
  for (int n = 0; n < mesh->num_nodes; ++n)
  {
    assert_int_equal(node_owners[node_offsets[n]], rank);
    for (int nn = node_offsets[n]+1; nn < node_offsets[n+1]; ++nn)
    {
      assert_true(node_owners[nn] >= 0);
      assert_true(node_owners[nn] < nprocs);
      assert_int_not_equal(node_owners[nn], rank);
    }
  }
  exchanger_fprintf(ex, stdout);
  exchanger_verify(ex, polymec_error);
  exchanger_free(ex);
  mesh_free(mesh);

  MPI_Barrier(MPI_COMM_WORLD);
}

void test_1v_node_exchanger_in_cube(void** state)
{
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  int nx = 8, ny = 8, nz = 8;
  mesh_t* mesh = create_uniform_mesh(MPI_COMM_SELF, nx, ny, nz, &bbox);
  partition_mesh(&mesh, MPI_COMM_WORLD, NULL, 0.0);
  exchanger_t* ex = mesh_1v_node_exchanger_new(mesh);
  exchanger_fprintf(ex, stdout);
  exchanger_verify(ex, polymec_error);
  exchanger_free(ex);
  mesh_free(mesh);

  MPI_Barrier(MPI_COMM_WORLD);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_nv_node_exchanger_on_line),
//    unit_test(test_1v_node_exchanger_on_line),
    unit_test(test_nv_node_exchanger_in_plane),
//    unit_test(test_1v_node_exchanger_in_plane),
    unit_test(test_nv_node_exchanger_in_cube),
//    unit_test(test_1v_node_exchanger_in_cube)
  };
  return run_tests(tests);
}
