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

// These tests are based on the solid body rotation test described in
// Nair and Jablonowski, "Moving Vortices on the Sphere: a Test Case for
// Horizontal Advection Problems", Mon. Wea. Rev. (2007).

// Horizontal solid body rotation velocity (u, v) as functions of lambda,
// theta.
//   u = u0 * (cos(theta) * cos(alpha) + sin(theta) * cos(lambda) * sin(alpha))
//   v = -u0 * sin(lambda) * sin(alpha)
// Above:
// * u and v are zonal and meridional wind velocities
// * lambda and theta are the longitude and latitude
// * alpha is the angle between the axis of rotation and the polar axis.
typedef struct
{
  real_t u0, alpha;
} sbr_t;

static void sbr_eval(void* context, point_t* x, real_t* f)
{
  sbr_t* sbr = context;
  real_t u0 = sbr->u0;
  real_t alpha = sbr->alpha;
  real_t lambda = x->x;
  real_t theta = x->y;
  // u = u0 * (cos(theta) * cos(alpha) + sin(theta) * cos(lambda) * sin(alpha))
  f[0] = u0*(cos(theta)*cos(alpha) + sin(theta)*cos(lambda)*sin(alpha));
  // v = -u0 * sin(lambda) * sin(alpha)
  f[1] = -u0*sin(lambda)*sin(alpha);
}

static sp_func_t* sbr_new(real_t u0, real_t alpha)
{
  sbr_t* sbr = polymec_malloc(sizeof(sbr_t));
  sbr->u0 = u0;
  sbr->alpha = alpha;
  sp_func_vtable vtable = {.eval = sbr_eval, .dtor = polymec_free};
  return sp_func_new("Solid body rotation", sbr, vtable,
                     SP_FUNC_HETEROGENEOUS, 2);
}

extern blockmesh_t* create_cubed_sphere(MPI_Comm comm,
                                        int block_nxy, int block_nz,
                                        int patch_nxy, int patch_nz,
                                        real_t R1, real_t R2);

static blockmesh_t* test_mesh(MPI_Comm comm)
{
  return create_cubed_sphere(MPI_COMM_SELF, 2, 2, 10, 10, 0.9, 1.0);
}

static void test_cell_field(void** state, blockmesh_t* mesh)
{
  sp_func_t* sbr = sbr_new(1.0, 0.0);
  blockmesh_field_t* f = blockmesh_field_new(mesh, UNIMESH_CELL, 2);
  blockmesh_field_apply(f, sbr);
  blockmesh_field_free(f);
  blockmesh_free(mesh);
  release_ref(sbr);
}

static void test_face_fields(void** state, blockmesh_t* mesh)
{
  sp_func_t* sbr = sbr_new(1.0, 0.0);
  blockmesh_field_t* fx = blockmesh_field_new(mesh, UNIMESH_XFACE, 2);
  blockmesh_field_t* fy = blockmesh_field_new(mesh, UNIMESH_YFACE, 2);
  blockmesh_field_t* fz = blockmesh_field_new(mesh, UNIMESH_ZFACE, 2);
  blockmesh_field_apply(fx, sbr);
  blockmesh_field_apply(fy, sbr);
  blockmesh_field_apply(fz, sbr);
  blockmesh_field_free(fx);
  blockmesh_field_free(fy);
  blockmesh_field_free(fz);
  blockmesh_free(mesh);
  release_ref(sbr);
}

static void test_edge_fields(void** state, blockmesh_t* mesh)
{
  sp_func_t* sbr = sbr_new(1.0, 0.0);
  blockmesh_field_t* fx = blockmesh_field_new(mesh, UNIMESH_XEDGE, 2);
  blockmesh_field_t* fy = blockmesh_field_new(mesh, UNIMESH_YEDGE, 2);
  blockmesh_field_t* fz = blockmesh_field_new(mesh, UNIMESH_ZEDGE, 2);
  blockmesh_field_apply(fx, sbr);
  blockmesh_field_apply(fy, sbr);
  blockmesh_field_apply(fz, sbr);
  blockmesh_field_free(fx);
  blockmesh_field_free(fy);
  blockmesh_field_free(fz);
  blockmesh_free(mesh);
}

static void test_node_field(void** state, blockmesh_t* mesh)
{
  sp_func_t* sbr = sbr_new(1.0, 0.0);
  blockmesh_field_t* f = blockmesh_field_new(mesh, UNIMESH_NODE, 2);
  blockmesh_field_apply(f, sbr);
  blockmesh_field_free(f);
  blockmesh_free(mesh);
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
