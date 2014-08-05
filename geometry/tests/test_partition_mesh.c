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
  real_t dx = 1.0/nx;
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = dx, .z1 = 0.0, .z2 = dx};
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
    exchanger_t* ex = mesh_exchanger(mesh);
    int pos = 0, proc, *indices, num_indices;
    int num_sends = 0, num_receives = 0;
    while (exchanger_next_send(ex, &pos, &proc, &indices, &num_indices))
    {
      for (int i = 0; i < num_indices; ++i)
        ++num_sends;
    }
    pos = 0;
    while (exchanger_next_receive(ex, &pos, &proc, &indices, &num_indices))
    {
      for (int i = 0; i < num_indices; ++i)
        ++num_receives;
    }
    assert_true((mesh->num_ghost_cells == 1) || (mesh->num_ghost_cells == 2));
    assert_true(num_sends == mesh->num_ghost_cells);
    assert_true(num_receives == mesh->num_ghost_cells);
  }
  else
    assert_int_equal(0, mesh->num_ghost_cells);

  // Check the geometry of the mesh.
  int cell_volumes_are_ok = 1;
  for (int c = 0; c < mesh->num_cells; ++c)
  {
    if (fabs(mesh->cell_volumes[c] - dx*dx*dx) > 1e-12)
    {
      cell_volumes_are_ok = 0;
      break; 
    }
  }
  int face_areas_are_ok = 1;
  for (int f = 0; f < mesh->num_faces; ++f)
  {
    if (fabs(mesh->face_areas[f] - dx*dx) > 1e-12)
    {
      face_areas_are_ok = 0;
      break; 
    }
  }
  MPI_Allreduce(&cell_volumes_are_ok, &cell_volumes_are_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&face_areas_are_ok, &cell_volumes_are_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  assert_true(cell_volumes_are_ok);
  assert_true(face_areas_are_ok);

  // Plot it.
  double p[mesh->num_cells];
  for (int c = 0; c < mesh->num_cells; ++c)
    p[c] = 1.0*rank;
  silo_file_t* silo = silo_file_new(mesh->comm, "linear_mesh_partition", "linear_mesh_partition", 1, 0, 0, 0.0);
  silo_file_write_mesh(silo, "mesh", mesh);
  silo_file_write_scalar_cell_field(silo, "rank", "mesh", p);
  silo_file_close(silo);

  // Clean up.
  mesh_free(mesh);

  // Superficially check that the file is okay.
  int num_files, num_procs;
  assert_true(silo_file_query("linear_mesh_partition", "linear_mesh_partition",
                              &num_files, &num_procs, NULL));
  assert_int_equal(1, num_files);
  assert_int_equal(nprocs, num_procs);
}

void test_partition_slab_mesh(void** state)
{
  // Create a 50x50x1 uniform mesh.
  int nx = 50, ny = 50, nz = 1;
  real_t dx = 1.0/nx;
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = dx};
  mesh_t* mesh = create_uniform_mesh(MPI_COMM_SELF, nx, ny, nz, &bbox);

  // Partition it.
  exchanger_t* distributor = partition_mesh(&mesh, MPI_COMM_WORLD, NULL, 0.05);
  exchanger_free(distributor);

  // Check the geometry of the mesh.
  int cell_volumes_are_ok = 1;
  for (int c = 0; c < mesh->num_cells; ++c)
  {
    if (fabs(mesh->cell_volumes[c] - dx*dx*dx) > 1e-12)
    {
      cell_volumes_are_ok = 0;
      break; 
    }
  }
  int face_areas_are_ok = 1;
  for (int f = 0; f < mesh->num_faces; ++f)
  {
    if (fabs(mesh->face_areas[f] - dx*dx) > 1e-12)
    {
      face_areas_are_ok = 0;
      break; 
    }
  }
  MPI_Allreduce(&cell_volumes_are_ok, &cell_volumes_are_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&face_areas_are_ok, &cell_volumes_are_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  assert_true(cell_volumes_are_ok);
  assert_true(face_areas_are_ok);

  // Plot it.
  int nprocs, rank;
  MPI_Comm_size(mesh->comm, &nprocs);
  MPI_Comm_rank(mesh->comm, &rank);
  double p[mesh->num_cells];
  for (int c = 0; c < mesh->num_cells; ++c)
    p[c] = 1.0*rank;
  silo_file_t* silo = silo_file_new(mesh->comm, "slab_mesh_partition", "slab_mesh_partition", 1, 0, 0, 0.0);
  silo_file_write_mesh(silo, "mesh", mesh);
  silo_file_write_scalar_cell_field(silo, "rank", "mesh", p);
  silo_file_close(silo);

  // Clean up.
  mesh_free(mesh);

  // Superficially check that the file is okay.
  int num_files, num_procs;
  assert_true(silo_file_query("slab_mesh_partition", "slab_mesh_partition",
                              &num_files, &num_procs, NULL));
  assert_int_equal(1, num_files);
  assert_int_equal(nprocs, num_procs);
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

  // Check the geometry of the mesh.
  real_t dx = 1.0/nx;
  int cell_volumes_are_ok = 1;
  for (int c = 0; c < mesh->num_cells; ++c)
  {
    if (fabs(mesh->cell_volumes[c] - dx*dx*dx) > 1e-12)
    {
      cell_volumes_are_ok = 0;
      break; 
    }
  }
  int face_areas_are_ok = 1;
  for (int f = 0; f < mesh->num_faces; ++f)
  {
    if (fabs(mesh->face_areas[f] - dx*dx) > 1e-12)
    {
      face_areas_are_ok = 0;
      break; 
    }
  }
  MPI_Allreduce(&cell_volumes_are_ok, &cell_volumes_are_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&face_areas_are_ok, &cell_volumes_are_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  assert_true(cell_volumes_are_ok);
  assert_true(face_areas_are_ok);

  // Plot it.
  int nprocs, rank;
  MPI_Comm_size(mesh->comm, &nprocs);
  MPI_Comm_rank(mesh->comm, &rank);
  double p[mesh->num_cells];
  for (int c = 0; c < mesh->num_cells; ++c)
    p[c] = 1.0*rank;
  silo_file_t* silo = silo_file_new(mesh->comm, "box_mesh_partition", "box_mesh_partition", 1, 0, 0, 0.0);
  silo_file_write_mesh(silo, "mesh", mesh);
  silo_file_write_scalar_cell_field(silo, "rank", "mesh", p);
  silo_file_close(silo);

  // Clean up.
  mesh_free(mesh);

  // Superficially check that the file is okay.
  int num_files, num_procs;
  assert_true(silo_file_query("box_mesh_partition", "box_mesh_partition",
                              &num_files, &num_procs, NULL));
  assert_int_equal(1, num_files);
  assert_int_equal(nprocs, num_procs);
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
