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
#include "core/mesh.h"
#include "geometry/create_uniform_mesh.h"

static void test_single_cell_mesh_no_topo(void** state)
{
  // Create a single hexahedron without topology.
  mesh_t* mesh = mesh_new_with_cell_type(MPI_COMM_SELF, 1, 0, 6, 8, 6, 4);
  assert_int_equal(1, mesh->num_cells);
  assert_int_equal(0, mesh->num_ghost_cells);
  assert_int_equal(6, mesh->num_faces);
  assert_int_equal(0, mesh->num_edges);
  assert_int_equal(8, mesh->num_nodes);
  mesh_free(mesh);
}

static void test_single_cell_mesh_serialization(void** state)
{
  // Create a single hexahedron.
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  mesh_t* mesh1 = create_uniform_mesh(MPI_COMM_SELF, 1, 1, 1, &bbox);

  // Now serialize!
  serializer_t* s = mesh_serializer();
  byte_array_t* buffer = byte_array_new();
  size_t offset = 0;
  serializer_write(s, mesh1, buffer, &offset);
  offset = 0;
  mesh_t* mesh2 = serializer_read(s, buffer, &offset);
  byte_array_free(buffer);

  assert_int_equal(1, mesh2->num_cells);
  assert_int_equal(0, mesh2->num_ghost_cells);
  assert_int_equal(6, mesh2->num_faces);
  assert_int_equal(8, mesh2->num_nodes);
  mesh_free(mesh2);
  mesh_free(mesh1);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_single_cell_mesh_no_topo),
    cmocka_unit_test(test_single_cell_mesh_serialization),
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
