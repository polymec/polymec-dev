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
#include "geometry/create_uniform_polymesh.h"
#include "geometry/crop_polymesh.h"
#include "geometry/cylinder_sd_func.h"
#include "geometry/plane_sd_func.h"
#include "geometry/intersection_sd_func.h"
#include "geometry/sphere_sd_func.h"

static void test_cylindrical_crop(void** state)
{
  // Create a cubic uniform mesh.
  int Nx = 10, Ny = 10, Nz = 10;
  real_t dz = 1.0/Nz;
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = -dz, .z2 = 1.0+dz};
  polymesh_t* mesh = create_uniform_polymesh(MPI_COMM_WORLD, Nx, Ny, Nz, &bbox);

  // Create a cropped mesh using a cylinder.
  polymesh_t* almost_cropped_mesh;
  {
    point_t O = {.x = 0.5, .y = 0.5, .z = 0.5};
    sd_func_t* cyl = cylinder_sd_func_new(&O, 0.5, INWARD_NORMAL);
    vector_t ntop = {.x = 0.0, .y = 0.0, .z = -1.0};
    point_t xtop = {.x = 0.0, .y = 0.0, .z = 1.0 + dz};
    vector_t nbot = {.x = 0.0, .y = 0.0, .z = 1.0};
    point_t xbot = {.x = 0.0, .y = 0.0, .z = 0.0 - dz};
    sd_func_t* ptop = plane_sd_func_new(&ntop, &xtop);
    sd_func_t* pbot = plane_sd_func_new(&nbot, &xbot);
    sd_func_t* surfaces[] = {cyl, ptop, pbot};
    sd_func_t* boundary = intersection_sd_func_new(surfaces, 3);
    almost_cropped_mesh = crop_polymesh(mesh, boundary, PROJECT_NODES);
    polymesh_free(mesh);
  }

  // Projecting the nodes *almost* works like it should, but there are some
  // artifacts. We chop off the top and bottom now to get rid of these.
  polymesh_t* cropped_mesh;
  {
    vector_t ntop = {.x = 0.0, .y = 0.0, .z = -1.0};
    point_t xtop = {.x = 0.0, .y = 0.0, .z = 1.0};
    vector_t nbot = {.x = 0.0, .y = 0.0, .z = 1.0};
    point_t xbot = {.x = 0.0, .y = 0.0, .z = 0.0};
    sd_func_t* ptop = plane_sd_func_new(&ntop, &xtop);
    sd_func_t* pbot = plane_sd_func_new(&nbot, &xbot);
    sd_func_t* surfaces[] = {ptop, pbot};
    sd_func_t* boundary = intersection_sd_func_new(surfaces, 2);
    cropped_mesh = crop_polymesh(almost_cropped_mesh, boundary, REMOVE_CELLS);
    polymesh_free(almost_cropped_mesh);
  }

  polymesh_free(cropped_mesh);
}

static void test_spherical_crop(void** state)
{
  // Create a cubic uniform mesh.
  int Nx = 10, Ny = 10, Nz = 10;
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  polymesh_t* mesh = create_uniform_polymesh(MPI_COMM_WORLD, Nx, Ny, Nz, &bbox);

  // Create a cropped mesh using a sphere.
  point_t O = {.x = 0.5, .y = 0.5, .z = 0.5};
  sd_func_t* boundary = sphere_sd_func_new(&O, 0.5, INWARD_NORMAL);
  polymesh_t* cropped_mesh = crop_polymesh(mesh, boundary, PROJECT_NODES);
  polymesh_free(mesh);
  polymesh_free(cropped_mesh);
}

int main(int argc, char* argv[])
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] =
  {
    cmocka_unit_test(test_cylindrical_crop),
    cmocka_unit_test(test_spherical_crop)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
