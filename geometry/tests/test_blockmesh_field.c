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
#include "core/options.h"
#include "geometry/blockmesh_field.h"

static blockmesh_t* test_mesh(MPI_Comm comm)
{
  return NULL;
}

static void test_cell_field(void** state, blockmesh_t* mesh)
{
}

static void test_face_fields(void** state, blockmesh_t* mesh)
{
}

static void test_edge_fields(void** state, blockmesh_t* mesh)
{
}

static void test_node_field(void** state, blockmesh_t* mesh)
{
}

static void test_serial_cell_field(void** state)
{
  blockmesh_t* mesh = test_mesh(MPI_COMM_SELF);
  test_cell_field(state, mesh);
}

static void test_serial_face_fields(void** state)
{
  blockmesh_t* mesh = test_mesh(MPI_COMM_SELF);
  test_face_fields(state, mesh);
}

static void test_serial_edge_fields(void** state)
{
  blockmesh_t* mesh = test_mesh(MPI_COMM_SELF);
  test_edge_fields(state, mesh);
}

static void test_serial_node_field(void** state)
{
  blockmesh_t* mesh = test_mesh(MPI_COMM_SELF);
  test_node_field(state, mesh);
}

static void test_parallel_cell_field(void** state)
{
  blockmesh_t* mesh = test_mesh(MPI_COMM_WORLD);
  test_cell_field(state, mesh);
}

static void test_parallel_face_fields(void** state)
{
  blockmesh_t* mesh = test_mesh(MPI_COMM_WORLD);
  test_face_fields(state, mesh);
}

static void test_parallel_edge_fields(void** state)
{
  blockmesh_t* mesh = test_mesh(MPI_COMM_WORLD);
  test_edge_fields(state, mesh);
}

static void test_parallel_node_field(void** state)
{
  blockmesh_t* mesh = test_mesh(MPI_COMM_WORLD);
  test_node_field(state, mesh);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_serial_cell_field),
    cmocka_unit_test(test_serial_face_fields),
    cmocka_unit_test(test_serial_edge_fields),
    cmocka_unit_test(test_serial_node_field),
    cmocka_unit_test(test_parallel_cell_field),
    cmocka_unit_test(test_parallel_face_fields),
    cmocka_unit_test(test_parallel_edge_fields),
    cmocka_unit_test(test_parallel_node_field)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
