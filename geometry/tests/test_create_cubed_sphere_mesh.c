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
#include "geometry/create_cubed_sphere_mesh.h"

void test_cubed_sphere_mesh(void** state, real_t r, real_t R)
{
  // Create a cubed sphere mesh.
  mesh_t* mesh = create_cubed_sphere_mesh(MPI_COMM_WORLD, 10, 10, r, R, "R");
  assert_true(mesh_verify_topology(mesh, polymec_error));
//  assert_int_equal(5000, mesh->num_cells);
  assert_true(mesh->comm == MPI_COMM_WORLD);

  char name[FILENAME_MAX];
  snprintf(name, FILENAME_MAX, "cubed_sphere_r=%g,R=%g", r, R);
  silo_file_t* silo = silo_file_new(MPI_COMM_WORLD, name, "", 1, 0, 0, 0.0);
  silo_file_write_mesh(silo, "mesh", mesh);
  silo_file_close(silo);

  // Clean up.
  mesh_free(mesh);
}

void test_create_cubed_sphere_mesh(void** state)
{
  test_cubed_sphere_mesh(state, 0.35, 0.5);
}

void test_cubed_spherical_shell_mesh(void** state, real_t r, real_t R)
{
  // Create a cubed spherical shell mesh.
  mesh_t* mesh = create_cubed_spherical_shell_mesh(MPI_COMM_WORLD, 10, 10, r, R,
                                                   "r1", "r2");
  assert_true(mesh_verify_topology(mesh, polymec_error));
//  assert_int_equal(8000, mesh->num_cells);

  char name[FILENAME_MAX];
  snprintf(name, FILENAME_MAX, "cubed_spherical_shell_r=%g,R=%g", r, R);
  silo_file_t* silo = silo_file_new(MPI_COMM_WORLD, name, "", 1, 0, 0, 0.0);
  silo_file_write_mesh(silo, "mesh", mesh);
  silo_file_close(silo);

  // Clean up.
  mesh_free(mesh);
}

void test_create_cubed_spherical_shell_mesh(void** state)
{
  test_cubed_spherical_shell_mesh(state, 0.35, 0.5);
}

void test_cubed_sphere_panel(void** state, real_t r, real_t R)
{
  // Create a cubed sphere panel.
  mesh_t* mesh = create_cubed_sphere_panel(MPI_COMM_WORLD, 10, 10, r, R, "north");
  assert_true(mesh_verify_topology(mesh, polymec_error));
//  assert_int_equal(8000, mesh->num_cells);

  char name[FILENAME_MAX];
  snprintf(name, FILENAME_MAX, "cubed_sphere_panel_north_r=%g,R=%g", r, R);
  silo_file_t* silo = silo_file_new(MPI_COMM_WORLD, name, "", 1, 0, 0, 0.0);
  silo_file_write_mesh(silo, "mesh", mesh);
  silo_file_close(silo);

  // Clean up.
  mesh_free(mesh);
}

void test_create_cubed_sphere_panel(void** state)
{
  test_cubed_sphere_panel(state, 0.35, 0.5);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
//    cmocka_unit_test(test_create_cubed_sphere_mesh),
    cmocka_unit_test(test_create_cubed_spherical_shell_mesh),
    cmocka_unit_test(test_create_cubed_sphere_panel)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
