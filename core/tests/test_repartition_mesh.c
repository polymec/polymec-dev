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
#include "core/repartition.h"
#include "geometry/create_uniform_mesh.h"

static void test_repartition_uniform_mesh_of_size(void** state, const char* prefix, int nx, int ny, int nz)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Create an nx x ny x nz uniform mesh.
  real_t dx = 1.0/MAX(MAX(1.0/nx, 1.0/ny), 1.0/nz);
  bbox_t bbox = {.x1 = 0.0, .x2 = nx*dx, .y1 = 0.0, .y2 = ny*dx, .z1 = 0.0, .z2 = nz*dx};
  mesh_t* mesh = create_uniform_mesh(MPI_COMM_WORLD, nx, ny, nz, &bbox);
  mesh_verify(mesh);

  // Repartition it.
  exchanger_t* migrator = repartition_mesh(&mesh, NULL, 0.05);
  exchanger_free(migrator);

  // Plot it.
  double r[mesh->num_cells];
  for (int c = 0; c < mesh->num_cells; ++c)
    r[c] = 1.0*rank;
  silo_file_t* silo = silo_file_new(mesh->comm, prefix, prefix, 1, 0, 0, 0.0);
  silo_file_write_mesh(silo, "mesh", mesh);
  silo_file_write_scalar_cell_field(silo, "rank", "mesh", r);
  silo_file_close(silo);

  // Clean up.
  mesh_free(mesh);
}

void test_repartition_x_uniform_mesh(void** state)
{
  test_repartition_uniform_mesh_of_size(state, "uniform_x_mesh_repartition", 128, 1, 1);
}

void test_repartition_xy_uniform_mesh(void** state)
{
  test_repartition_uniform_mesh_of_size(state, "uniform_xy_mesh_repartition", 128, 128, 1);
}

void test_repartition_xyz_uniform_mesh(void** state)
{
  test_repartition_uniform_mesh_of_size(state, "uniform_xyz_mesh_repartition", 32, 32, 32);
}

void test_repartition_tiny_x_uniform_mesh(void** state)
{
  test_repartition_uniform_mesh_of_size(state, "tiny_x_mesh_repartition", 4, 1, 1);
}

void test_repartition_tiny_xy_uniform_mesh(void** state)
{
  test_repartition_uniform_mesh_of_size(state, "tiny_xy_mesh_repartition", 4, 4, 1);
}

void test_repartition_tiny_xyz_uniform_mesh(void** state)
{
  test_repartition_uniform_mesh_of_size(state, "tiny_xyz_mesh_repartition", 4, 4, 4);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_repartition_tiny_x_uniform_mesh),
    unit_test(test_repartition_tiny_xy_uniform_mesh),
//    unit_test(test_repartition_tiny_xyz_uniform_mesh),
    unit_test(test_repartition_x_uniform_mesh),
    unit_test(test_repartition_xy_uniform_mesh)
//    unit_test(test_repartition_xyz_uniform_mesh)
  };
  return run_tests(tests);
}
