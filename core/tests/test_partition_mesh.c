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
#include "core/silo_file.h"
#include "core/partition_mesh.h"
#include "geometry/create_uniform_mesh.h"

#if POLYMEC_HAVE_DOUBLE_PRECISION
static const real_t volume_tolerance = 1e-12;
static const real_t area_tolerance = 1e-12;
#else
static const real_t volume_tolerance = 1e-6;
static const real_t area_tolerance = 1e-6;
#endif
static void test_partition_linear_mesh(void** state)
{
  // Create a 100x1x1 uniform mesh.
  int nx = 100, ny = 1, nz = 1;
  real_t dx = 1.0/nx;
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = dx, .z1 = 0.0, .z2 = dx};
  mesh_t* mesh = create_uniform_mesh(MPI_COMM_SELF, nx, ny, nz, &bbox);

  // Partition it.
  migrator_t* m = partition_mesh(&mesh, MPI_COMM_WORLD, NULL, 0.05);
  migrator_verify(m, polymec_error);
  m = NULL;

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
      num_sends += num_indices;
    pos = 0;
    while (exchanger_next_receive(ex, &pos, &proc, &indices, &num_indices))
      num_receives += num_indices;
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
    if (!reals_nearly_equal(mesh->cell_volumes[c], dx*dx*dx, volume_tolerance))
    {
      cell_volumes_are_ok = 0;
      break; 
    }
  }
  int face_areas_are_ok = 1;
  for (int f = 0; f < mesh->num_faces; ++f)
  {
    if (!reals_nearly_equal(mesh->face_areas[f], dx*dx, area_tolerance))
    {
      face_areas_are_ok = 0;
      break; 
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, &cell_volumes_are_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &face_areas_are_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  assert_true(cell_volumes_are_ok);
  assert_true(face_areas_are_ok);

  // Check the resulting exchanger.
  exchanger_verify(mesh_exchanger(mesh), polymec_error);

  // Plot it.
  real_t p[mesh->num_cells];
  for (int c = 0; c < mesh->num_cells; ++c)
    p[c] = 1.0*rank;
  silo_file_t* silo = silo_file_new(mesh->comm, "linear_mesh_partition", "linear_mesh_partition", 1, 0, 0, 0.0);
  silo_file_write_mesh(silo, "mesh", mesh);
  silo_field_metadata_t* p_metadata = silo_field_metadata_new();
  silo_field_metadata_set_label(p_metadata, "P");
  silo_field_metadata_set_conserved(p_metadata, false);
  silo_file_write_scalar_cell_field(silo, "rank", "mesh", p, p_metadata);
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

static void test_partition_slab_mesh(void** state)
{
  // Create a 50x50x1 uniform mesh.
  int nx = 50, ny = 50, nz = 1;
  real_t dx = 1.0/nx;
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = dx};
  mesh_t* mesh = create_uniform_mesh(MPI_COMM_SELF, nx, ny, nz, &bbox);

  // Partition it.
  migrator_t* m = partition_mesh(&mesh, MPI_COMM_WORLD, NULL, 0.05);
  m = NULL;

  // Check the geometry of the mesh.
  int cell_volumes_are_ok = 1;
  for (int c = 0; c < mesh->num_cells; ++c)
  {
    if (!reals_nearly_equal(mesh->cell_volumes[c], dx*dx*dx, volume_tolerance))
    {
      cell_volumes_are_ok = 0;
      break; 
    }
  }
  int face_areas_are_ok = 1;
  for (int f = 0; f < mesh->num_faces; ++f)
  {
    if (!reals_nearly_equal(mesh->face_areas[f], dx*dx, area_tolerance))
    {
      face_areas_are_ok = 0;
      break; 
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, &cell_volumes_are_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &face_areas_are_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  assert_true(cell_volumes_are_ok);
  assert_true(face_areas_are_ok);

  // Plot it.
  int nprocs, rank;
  MPI_Comm_size(mesh->comm, &nprocs);
  MPI_Comm_rank(mesh->comm, &rank);
  real_t p[mesh->num_cells];
  for (int c = 0; c < mesh->num_cells; ++c)
    p[c] = 1.0*rank;
  silo_file_t* silo = silo_file_new(mesh->comm, "slab_mesh_partition", "slab_mesh_partition", 1, 0, 0, 0.0);
  silo_file_write_mesh(silo, "mesh", mesh);
  silo_file_write_scalar_cell_field(silo, "rank", "mesh", p, NULL);
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

static void test_partition_box_mesh(void** state)
{
  // Create a 20x20x20 uniform mesh.
  int nx = 20, ny = 20, nz = 20;
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  mesh_t* mesh = create_uniform_mesh(MPI_COMM_SELF, nx, ny, nz, &bbox);

  // Partition it.
  migrator_t* m = partition_mesh(&mesh, MPI_COMM_WORLD, NULL, 0.05);
  m = NULL;

  // Check the geometry of the mesh.
  real_t dx = 1.0/nx;
  int cell_volumes_are_ok = 1;
  for (int c = 0; c < mesh->num_cells; ++c)
  {
    if (!reals_nearly_equal(mesh->cell_volumes[c], dx*dx*dx, volume_tolerance))
    {
      cell_volumes_are_ok = 0;
      break; 
    }
  }
  int face_areas_are_ok = 1;
  for (int f = 0; f < mesh->num_faces; ++f)
  {
    if (!reals_nearly_equal(mesh->face_areas[f], dx*dx, area_tolerance))
    {
      face_areas_are_ok = 0;
      break; 
    }
  }
  MPI_Allreduce(MPI_IN_PLACE, &cell_volumes_are_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &face_areas_are_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
  assert_true(cell_volumes_are_ok);
  assert_true(face_areas_are_ok);

  // Plot it.
  int nprocs, rank;
  MPI_Comm_size(mesh->comm, &nprocs);
  MPI_Comm_rank(mesh->comm, &rank);
  real_t p[mesh->num_cells];
  for (int c = 0; c < mesh->num_cells; ++c)
    p[c] = 1.0*rank;
  silo_file_t* silo = silo_file_new(mesh->comm, "box_mesh_partition", "box_mesh_partition", 1, 0, 0, 0.0);
  silo_file_write_mesh(silo, "mesh", mesh);
  silo_file_write_scalar_cell_field(silo, "rank", "mesh", p, NULL);
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
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_partition_linear_mesh),
    cmocka_unit_test(test_partition_slab_mesh),
    cmocka_unit_test(test_partition_box_mesh)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
