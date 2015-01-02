// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this 
// list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice, 
// this list of conditions and the following disclaimer in the documentation 
// and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmockery.h"
#include "core/silo_file.h"
#include "geometry/create_uniform_mesh.h"
#include "geometry/crop_mesh.h"
#include "geometry/cylinder.h"
#include "geometry/plane.h"
#include "geometry/intersection.h"
#include "geometry/sphere.h"

void test_cylindrical_crop(void** state)
{
  // Create a cubic uniform mesh.
  int Nx = 10, Ny = 10, Nz = 10;
  real_t dz = 1.0/Nz;
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = -dz, .z2 = 1.0+dz};
  mesh_t* mesh = create_uniform_mesh(MPI_COMM_WORLD, Nx, Ny, Nz, &bbox);

  // Create a cropped mesh using a cylinder.
  mesh_t* almost_cropped_mesh;
  {
    point_t O = {.x = 0.5, .y = 0.5, .z = 0.5};
    sp_func_t* cyl = cylinder_new(&O, 0.5, INWARD_NORMAL);
    vector_t ntop = {.x = 0.0, .y = 0.0, .z = -1.0}; 
    point_t xtop = {.x = 0.0, .y = 0.0, .z = 1.0 + dz};
    vector_t nbot = {.x = 0.0, .y = 0.0, .z = 1.0}; 
    point_t xbot = {.x = 0.0, .y = 0.0, .z = 0.0 - dz};
    sp_func_t* ptop = plane_new(&ntop, &xtop);
    sp_func_t* pbot = plane_new(&nbot, &xbot);
    sp_func_t* surfaces[] = {cyl, ptop, pbot};
    sp_func_t* boundary = intersection_new(surfaces, 3);
    almost_cropped_mesh = crop_mesh(mesh, boundary, PROJECT_NODES);
    mesh_free(mesh);
  }

  // Projecting the nodes *almost* works like it should, but there are some 
  // artifacts. We chop off the top and bottom now to get rid of these.
  mesh_t* cropped_mesh;
  {
    vector_t ntop = {.x = 0.0, .y = 0.0, .z = -1.0}; 
    point_t xtop = {.x = 0.0, .y = 0.0, .z = 1.0};
    vector_t nbot = {.x = 0.0, .y = 0.0, .z = 1.0}; 
    point_t xbot = {.x = 0.0, .y = 0.0, .z = 0.0};
    sp_func_t* ptop = plane_new(&ntop, &xtop);
    sp_func_t* pbot = plane_new(&nbot, &xbot);
    sp_func_t* surfaces[] = {ptop, pbot};
    sp_func_t* boundary = intersection_new(surfaces, 2);
    cropped_mesh = crop_mesh(almost_cropped_mesh, boundary, REMOVE_CELLS);
    mesh_free(almost_cropped_mesh);
  }

  // Plot the cropped mesh.
  double ones[Nx*Ny*Nz];
  for (int c = 0; c < Nx*Ny*Nz; ++c)
    ones[c] = 1.0*c;
  silo_file_t* silo = silo_file_new(cropped_mesh->comm, "cyl_cropped_mesh", "", 1, 0, 0, 0.0);
  silo_file_write_mesh(silo, "mesh", cropped_mesh);
  silo_file_write_scalar_cell_field(silo, "solution", "mesh", ones);
  silo_file_close(silo);

  mesh_free(cropped_mesh);
}

void test_spherical_crop(void** state)
{
  // Create a cubic uniform mesh.
  int Nx = 10, Ny = 10, Nz = 10;
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  mesh_t* mesh = create_uniform_mesh(MPI_COMM_WORLD, Nx, Ny, Nz, &bbox);

  // Create a cropped mesh using a sphere.
  point_t O = {.x = 0.5, .y = 0.5, .z = 0.5};
  sp_func_t* boundary = sphere_new(&O, 0.5, INWARD_NORMAL);
  mesh_t* cropped_mesh = crop_mesh(mesh, boundary, PROJECT_NODES);
  mesh_free(mesh);

  // Plot the cropped mesh.
  double ones[Nx*Ny*Nz];
  for (int c = 0; c < Nx*Ny*Nz; ++c)
    ones[c] = 1.0*c;
  silo_file_t* silo = silo_file_new(cropped_mesh->comm, "sph_cropped_mesh", "", 1, 0, 0, 0.0);
  silo_file_write_mesh(silo, "mesh", cropped_mesh);
  silo_file_write_scalar_cell_field(silo, "solution", "mesh", ones);
  silo_file_close(silo);

  mesh_free(cropped_mesh);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_cylindrical_crop),
    unit_test(test_spherical_crop)
  };
  return run_tests(tests);
}
