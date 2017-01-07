// Copyright (c) 2012-2017, Jeffrey N. Johnson
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
#include "geometry/create_cubed_cylinder_mesh.h"

void test_cubed_cylinder_mesh(void** state, real_t R, real_t L, real_t l, bool curved)
{
  // Create a cubed cylinder mesh with a square center block.
  mesh_t* mesh = create_cubed_cylinder_mesh(MPI_COMM_WORLD, 10, 20, R, L, l, curved,
                                            "R", "bottom", "top");
  assert_true(mesh_verify_topology(mesh, polymec_error));
//  assert_int_equal(5000, mesh->num_cells);
  assert_true(mesh->comm == MPI_COMM_WORLD);

  char name[FILENAME_MAX];
  if (curved)
    snprintf(name, FILENAME_MAX, "cubed_circular_cylinder_R=%g,L=%g,l=%g", R, L, l);
  else
    snprintf(name, FILENAME_MAX, "cubed_cylinder_R=%g,L=%g,l=%g", R, L, l);
  silo_file_t* silo = silo_file_new(MPI_COMM_WORLD, name, "", 1, 0, 0, 0.0);
  silo_file_write_mesh(silo, "mesh", mesh);
  silo_file_close(silo);

  // Clean up.
  mesh_free(mesh);
}

void test_create_cubed_cylinder_mesh(void** state)
{
  test_cubed_cylinder_mesh(state, 0.5, 1.0, 0.35, false);
}

void test_create_circular_cubed_cylinder_mesh(void** state)
{
  test_cubed_cylinder_mesh(state, 0.5, 1.0, 0.35, true);
}

void test_cubed_cylindrical_shell_mesh(void** state, real_t r, real_t R, real_t L)
{
  // Create a cubed cylindrcal shell mesh.
  mesh_t* mesh = create_cubed_cylindrical_shell_mesh(MPI_COMM_WORLD, 10, 20, r, R, L,
                                                     "r1", "r2", "bottom", "top");
  assert_true(mesh_verify_topology(mesh, polymec_error));
  assert_int_equal(8000, mesh->num_cells);

  char name[FILENAME_MAX];
  snprintf(name, FILENAME_MAX, "cubed_cylindrical_shell_r=%g,R=%g,L=%g", r, R, L);
  silo_file_t* silo = silo_file_new(MPI_COMM_WORLD, name, "", 1, 0, 0, 0.0);
  silo_file_write_mesh(silo, "mesh", mesh);
  silo_file_close(silo);

  // Clean up.
  mesh_free(mesh);
}

void test_create_cubed_cylindrical_shell_mesh(void** state)
{
  test_cubed_cylindrical_shell_mesh(state, 0.35, 0.5, 1.0);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_create_cubed_cylinder_mesh),
    cmocka_unit_test(test_create_circular_cubed_cylinder_mesh),
    cmocka_unit_test(test_create_cubed_cylindrical_shell_mesh)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
