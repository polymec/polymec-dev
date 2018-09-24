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
#include "geometry/create_uniform_polymesh.h"
#include "io/silo_file.h"

static void test_plot_uniform_mesh_with_num_files(void** state, int num_files)
{
  // Create a uniform mesh.
  int nx = 10, ny = 10, nz = 10;
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  polymesh_t* mesh = create_uniform_polymesh(MPI_COMM_WORLD, nx, ny, nz, &bbox);

  // Plot it and some field data.
  char prefix[FILENAME_MAX];
  snprintf(prefix, FILENAME_MAX, "uniform_mesh_%dx%dx%d_%df", nx, ny, nz, num_files);
  silo_file_t* silo = silo_file_new(mesh->comm, prefix, "", num_files, 0, 0.0);
  silo_file_write_polymesh(silo, "mesh", mesh);

  polymesh_field_t* cfield = polymesh_field_new(mesh, POLYMESH_CELL, 1);
  DECLARE_POLYMESH_FIELD_ARRAY(cvals, cfield);
  const char* cnames[] = {"solution"};
  for (int c = 0; c < 4*4*4; ++c)
    cvals[c][0] = 1.0*c;
  silo_file_write_polymesh_field(silo, cnames, "mesh", cfield, NULL);
  silo_file_close(silo);

  // Query the plot file to make sure its numbers are good.
  int nprocs, my_num_files, num_mpi_procs;
  MPI_Comm_size(mesh->comm, &nprocs);
  char dir_name[FILENAME_MAX];
  if (nprocs > 1)
    snprintf(dir_name, FILENAME_MAX, "%s_%dprocs", prefix, nprocs);
  else
    snprintf(dir_name, FILENAME_MAX, ".");
  assert_true(silo_file_query(prefix, dir_name, &my_num_files, &num_mpi_procs, NULL));
  assert_int_equal(num_files, my_num_files);
  assert_int_equal(nprocs, num_mpi_procs);

  // Get cycles too.
  int_slist_t* cycles = int_slist_new();
  assert_true(silo_file_query(prefix, dir_name, &my_num_files, &num_mpi_procs, cycles));
  assert_int_equal(num_files, my_num_files);
  assert_int_equal(nprocs, num_mpi_procs);
  assert_int_equal(1, cycles->size);
  if (cycles->size > 0)
    assert_int_equal(0, cycles->front->value);
  int_slist_free(cycles);

  // Clean up.
  polymesh_free(mesh);
}

static void test_plot_uniform_mesh_to_single_file(void** state)
{
  test_plot_uniform_mesh_with_num_files(state, 1);
}

static void test_plot_uniform_mesh_to_n_files(void** state)
{
  int nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  test_plot_uniform_mesh_with_num_files(state, nprocs);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_plot_uniform_mesh_to_single_file),
    cmocka_unit_test(test_plot_uniform_mesh_to_n_files)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
