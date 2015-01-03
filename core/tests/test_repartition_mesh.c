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
#include "core/silo_file.h"
#include "core/partition_mesh.h"
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

  // Since the mesh is uniform, we can check the properties of each cell.
  for (int c = 0; c < mesh->num_cells; ++c)
  {
    assert_int_equal(6, mesh_cell_num_faces(mesh, c));
    real_t V = dx * dx * dx;
printf("%d: V[%d] = %g, should be %g\n", rank, c, mesh->cell_volumes[c], V);
    assert_true(fabs(mesh->cell_volumes[c] - V)/V < 1e-14);
  }

  // We can also check the properties of each face.
  for (int f = 0; f < mesh->num_faces; ++f)
  {
    assert_int_equal(4, mesh_face_num_edges(mesh, f));
    assert_int_equal(4, mesh_face_num_nodes(mesh, f));
    real_t A = dx * dx;
    assert_true(fabs(mesh->face_areas[f] - A)/A < 1e-14);
  }

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
