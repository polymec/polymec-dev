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

void test_create_uniform_mesh(void** state)
{
  // Create a 10x10x10 uniform mesh.
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  mesh_t* mesh = create_uniform_mesh(MPI_COMM_WORLD, 10, 10, 10, &bbox);
  mesh_verify(mesh);

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

void test_plot_uniform_mesh_with_num_files(void** state, int num_files)
{
  // Create a uniform mesh.
  int nx = 10, ny = 10, nz = 10;
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  mesh_t* mesh = create_uniform_mesh(MPI_COMM_WORLD, nx, ny, nz, &bbox);

  // Plot it.
  double ones[nx*ny*nz];
  for (int c = 0; c < ny*ny*nz; ++c)
    ones[c] = 1.0*c;
  char prefix[FILENAME_MAX];
  snprintf(prefix, FILENAME_MAX, "uniform_mesh_%dx%dx%d", nx, ny, nz);
  silo_file_t* silo = silo_file_new(mesh->comm, prefix, "", num_files, 0, 0, 0.0);
  silo_file_write_mesh(silo, "mesh", mesh);
  silo_file_write_scalar_cell_field(silo, "solution", "mesh", ones);
  silo_file_close(silo);

  // Query the plot file to make sure its numbers are good.
  int nprocs, my_num_files, num_mpi_procs;
  MPI_Comm_size(mesh->comm, &nprocs);
  char dir_name[FILENAME_MAX];
  if (nprocs > 1)
    snprintf(dir_name, FILENAME_MAX, "%s_%dprocs", prefix, nprocs);
  else
    snprintf(dir_name, FILENAME_MAX, ".");
  silo_file_query(prefix, dir_name, &my_num_files, &num_mpi_procs, NULL);
  assert_int_equal(num_files, my_num_files);
  assert_int_equal(nprocs, num_mpi_procs);

  // Get cycles too.
  int_slist_t* cycles = int_slist_new();
  silo_file_query(prefix, dir_name, &my_num_files, &num_mpi_procs, cycles);
  assert_int_equal(num_files, my_num_files);
  assert_int_equal(nprocs, num_mpi_procs);
  assert_int_equal(1, cycles->size);
  if (cycles->size > 0)
    assert_int_equal(0, cycles->front->value);
  int_slist_free(cycles);

  // Clean up.
  mesh_free(mesh);
}

void test_plot_uniform_mesh_to_single_file(void** state)
{
  test_plot_uniform_mesh_with_num_files(state, 1);
}

void test_plot_uniform_mesh_to_n_files(void** state)
{
  int nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  test_plot_uniform_mesh_with_num_files(state, nprocs);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_create_uniform_mesh),
    unit_test(test_plot_uniform_mesh_to_single_file),
    unit_test(test_plot_uniform_mesh_to_n_files)
  };
  return run_tests(tests);
}
