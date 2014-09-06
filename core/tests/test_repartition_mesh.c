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

static void test_repartition_uniform_mesh_of_size(void** state, int nx, int ny, int nz)
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Create an nx x ny x nz uniform mesh.
  real_t dx = 1.0/MAX(MAX(1.0/nx, 1.0/ny), 1.0/nz);
  bbox_t bbox = {.x1 = 0.0, .x2 = nx*dx, .y1 = 0.0, .y2 = ny*dx, .z1 = 0.0, .z2 = nz*dx};
  mesh_t* mesh = create_uniform_mesh(MPI_COMM_WORLD, nx, ny, nz, &bbox);
  mesh_verify_topology(mesh, polymec_error);

  // Repartition it.
  exchanger_t* migrator = repartition_mesh(&mesh, NULL, 0.05);
  exchanger_free(migrator);

  // Check the total volume.
  real_t V = 0.0, V_actual = nx*dx * ny*dx * ny*dx;
  for (int c = 0; c < mesh->num_cells; ++c)
    V += mesh->cell_volumes[c];
  real_t V_total;
  MPI_Allreduce(&V, &V_total, 1, MPI_REAL_T, MPI_SUM, MPI_COMM_WORLD);
printf("%d: V_total = %g, V_actual = %g, |V_total - V_actual| = %g\n", rank, V_total, V_actual, fabs(V_total - V_actual));
  assert_true(fabs(V_total - V_actual) < 1e-14);

  // Plot it.
  double r[mesh->num_cells];
  for (int c = 0; c < mesh->num_cells; ++c)
    r[c] = 1.0*rank;
  char prefix[FILENAME_MAX];
  snprintf(prefix, FILENAME_MAX, "%dx%dx%d_uniform_mesh_repartition", nx, ny, nz);
  silo_file_t* silo = silo_file_new(mesh->comm, prefix, prefix, 1, 0, 0, 0.0);
  silo_file_write_mesh(silo, "mesh", mesh);
  silo_file_write_scalar_cell_field(silo, "rank", "mesh", r);
  silo_file_close(silo);

  // Clean up.
  mesh_free(mesh);
}

void test_repartition_4x1x1_uniform_mesh(void** state)
{
  test_repartition_uniform_mesh_of_size(state, 4, 1, 1);
}

void test_repartition_2x2x1_uniform_mesh(void** state)
{
  test_repartition_uniform_mesh_of_size(state, 2, 2, 1);
}

void test_repartition_4x4x1_uniform_mesh(void** state)
{
  test_repartition_uniform_mesh_of_size(state, 4, 4, 1);
}

void test_repartition_2x2x2_uniform_mesh(void** state)
{
  test_repartition_uniform_mesh_of_size(state, 2, 2, 2);
}

void test_repartition_4x4x4_uniform_mesh(void** state)
{
  test_repartition_uniform_mesh_of_size(state, 4, 4, 4);
}

void test_repartition_128x1x1_uniform_mesh(void** state)
{
  test_repartition_uniform_mesh_of_size(state, 128, 1, 1);
}

void test_repartition_128x128x1_uniform_mesh(void** state)
{
  test_repartition_uniform_mesh_of_size(state, 128, 128, 1);
}

void test_repartition_32x32x32_uniform_mesh(void** state)
{
  test_repartition_uniform_mesh_of_size(state, 32, 32, 32);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_repartition_4x1x1_uniform_mesh),
    unit_test(test_repartition_2x2x1_uniform_mesh),
    unit_test(test_repartition_4x4x1_uniform_mesh),
    unit_test(test_repartition_2x2x2_uniform_mesh),
    unit_test(test_repartition_4x4x4_uniform_mesh),
    unit_test(test_repartition_128x1x1_uniform_mesh),
    unit_test(test_repartition_128x128x1_uniform_mesh),
    unit_test(test_repartition_32x32x32_uniform_mesh)
  };
  return run_tests(tests);
}
