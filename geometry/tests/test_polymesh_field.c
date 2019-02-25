// Copyright (c) 2012-2019, Jeffrey N. Johnson
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
#include "geometry/polymesh_field.h"
#include "geometry/create_uniform_polymesh.h"

static void test_single_cell_mesh_no_topo(void** state)
{
  // Create a single hexahedron without topology.
  polymesh_t* mesh = polymesh_new_with_cell_type(MPI_COMM_SELF, 1, 0, 6, 8, 6, 4);
  polymesh_construct_edges(mesh);
  for (int c = 1; c < 4; ++c)
  {
    polymesh_field_t* fc = polymesh_field_new(mesh, POLYMESH_CELL, c);
    assert_true(fc->centering == POLYMESH_CELL);
    assert_int_equal(c, fc->num_components);
    assert_int_equal(mesh->num_cells, fc->num_local_values);
    assert_int_equal(mesh->num_ghost_cells, fc->num_ghost_values);
    polymesh_field_set_exchanger(fc, polymesh_exchanger(mesh, POLYMESH_CELL));
    polymesh_field_free(fc);

    polymesh_field_t* ff = polymesh_field_new(mesh, POLYMESH_FACE, c);
    assert_true(ff->centering == POLYMESH_FACE);
    assert_int_equal(c, ff->num_components);
    assert_int_equal(mesh->num_faces, ff->num_local_values);
    assert_int_equal(0, ff->num_ghost_values);
    polymesh_field_set_exchanger(ff, polymesh_exchanger(mesh, POLYMESH_FACE));
    polymesh_field_free(ff);

    polymesh_field_t* fe = polymesh_field_new(mesh, POLYMESH_EDGE, c);
    assert_true(fe->centering == POLYMESH_EDGE);
    assert_int_equal(c, fe->num_components);
    assert_int_equal(mesh->num_edges, fe->num_local_values);
    assert_int_equal(0, fe->num_ghost_values);
    polymesh_field_set_exchanger(fe, polymesh_exchanger(mesh, POLYMESH_EDGE));
    polymesh_field_free(fe);

    polymesh_field_t* fn = polymesh_field_new(mesh, POLYMESH_NODE, c);
    assert_true(fn->centering == POLYMESH_NODE);
    assert_int_equal(c, fn->num_components);
    assert_int_equal(mesh->num_nodes, fn->num_local_values);
    assert_int_equal(0, fn->num_ghost_values);
    polymesh_field_set_exchanger(fn, polymesh_exchanger(mesh, POLYMESH_NODE));
    polymesh_field_free(fn);
  }
  polymesh_free(mesh);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_single_cell_mesh_no_topo)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
