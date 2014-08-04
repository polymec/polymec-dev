// Copyright (c) 2012-2014, Jeffrey N. Johnson
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this 
// list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice, 
// this list of conditions and the following disclaimer in the documentation 
// and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmockery.h"
#include "core/silo_file.h"
#include "geometry/create_uniform_mesh.h"
#include "geometry/repartition.h"

void test_partition_linear_mesh(void** state)
{
  // Create a 100x1x1 uniform mesh.
  int nx = 100, ny = 1, nz = 1;
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  mesh_t* mesh = create_uniform_mesh(MPI_COMM_SELF, nx, ny, nz, &bbox);

  // Partition it.
  exchanger_t* distributor = partition_mesh(&mesh, MPI_COMM_WORLD, NULL, 0.05);
  exchanger_free(distributor);

  // Check the ghost cells.
  int rank, nprocs;
  MPI_Comm_rank(mesh->comm, &rank);
  MPI_Comm_size(mesh->comm, &nprocs);
  if (nprocs > 1)
  {
    if ((rank == 0) || (rank == (nprocs-1)))
      assert_int_equal(1, mesh->num_ghost_cells);
    else
      assert_int_equal(2, mesh->num_ghost_cells);
  }
  else
    assert_int_equal(0, mesh->num_ghost_cells);

  // Plot it.
  double ones[mesh->num_cells];
  for (int c = 0; c < mesh->num_cells; ++c)
    ones[c] = 1.0*rank;
  silo_file_t* silo = silo_file_new(mesh->comm, "linear_mesh_partition", "linear_mesh_partition", 1, 0, 0, 0.0);
  silo_file_write_mesh(silo, "mesh", mesh);
  silo_file_write_scalar_cell_field(silo, "rank", "mesh", ones);
  silo_file_close(silo);

  // Clean up.
  mesh_free(mesh);
}

void test_partition_slab_mesh(void** state)
{
  // Create a 50x50x1 uniform mesh.
  int nx = 50, ny = 50, nz = 1;
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  mesh_t* mesh = create_uniform_mesh(MPI_COMM_SELF, nx, ny, nz, &bbox);

  // Partition it.
  exchanger_t* distributor = partition_mesh(&mesh, MPI_COMM_WORLD, NULL, 0.05);
  exchanger_free(distributor);

  // Plot it.
  int rank;
  MPI_Comm_rank(mesh->comm, &rank);
  double ones[mesh->num_cells];
  for (int c = 0; c < mesh->num_cells; ++c)
    ones[c] = 1.0*rank;
  silo_file_t* silo = silo_file_new(mesh->comm, "slab_mesh_partition", "slab_mesh_partition", 1, 0, 0, 0.0);
  silo_file_write_mesh(silo, "mesh", mesh);
  silo_file_write_scalar_cell_field(silo, "rank", "mesh", ones);
  silo_file_close(silo);

  // Clean up.
  mesh_free(mesh);
}

void test_partition_box_mesh(void** state)
{
  // Create a 20x20x20 uniform mesh.
  int nx = 20, ny = 20, nz = 20;
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  mesh_t* mesh = create_uniform_mesh(MPI_COMM_SELF, nx, ny, nz, &bbox);

  // Partition it.
  exchanger_t* distributor = partition_mesh(&mesh, MPI_COMM_WORLD, NULL, 0.05);
  exchanger_free(distributor);

  // Plot it.
  int rank;
  MPI_Comm_rank(mesh->comm, &rank);
  double ones[mesh->num_cells];
  for (int c = 0; c < mesh->num_cells; ++c)
    ones[c] = 1.0*rank;
  silo_file_t* silo = silo_file_new(mesh->comm, "box_mesh_partition", "box_mesh_partition", 1, 0, 0, 0.0);
  silo_file_write_mesh(silo, "mesh", mesh);
  silo_file_write_scalar_cell_field(silo, "rank", "mesh", ones);
  silo_file_close(silo);

  // Clean up.
  mesh_free(mesh);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_partition_linear_mesh),
    unit_test(test_partition_slab_mesh),
    unit_test(test_partition_box_mesh)
  };
  return run_tests(tests);
}
