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
#include "geometry/create_cubed_cylinder_mesh.h"

void test_create_cubed_cylinder_mesh(void** state)
{
  // Create a cubed cylinder mesh with a square center block.
  real_t R = 0.5, L = 1.0;
  real_t l = 0.35, k = 0.0;
  mesh_t* mesh = create_cubed_cylinder_mesh(MPI_COMM_WORLD, 10, 10, R, L, l, k);
  assert_true(mesh_verify_topology(mesh, polymec_error));
//  assert_int_equal(5000, mesh->num_cells);
  assert_true(mesh->comm == MPI_COMM_WORLD);

  silo_file_t* silo = silo_file_new(MPI_COMM_WORLD, "cubed_cylinder", "", 1, 0, 0, 0.0);
  silo_file_write_mesh(silo, "mesh", mesh);
  silo_file_close(silo);

  // Clean up.
  mesh_free(mesh);
}

void test_create_cubed_cylindrical_shell_mesh(void** state)
{
  // Create a cubed cylindrcal shell mesh.
  real_t r = 0.25, R = 0.5, L = 1.0;
  mesh_t* mesh = create_cubed_cylindrical_shell_mesh(MPI_COMM_WORLD, 10, 10, r, R, L);
  assert_true(mesh_verify_topology(mesh, polymec_error));
  assert_int_equal(4000, mesh->num_cells);

  silo_file_t* silo = silo_file_new(MPI_COMM_WORLD, "cubed_cylindrical_shell", "", 1, 0, 0, 0.0);
  silo_file_write_mesh(silo, "mesh", mesh);
  silo_file_close(silo);

  // Clean up.
  mesh_free(mesh);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_create_cubed_cylinder_mesh)
//    unit_test(test_create_cubed_cylindrical_shell_mesh)
  };
  return run_tests(tests);
}
