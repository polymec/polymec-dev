// Copyright (c) 2012-2013, Jeffrey N. Johnson
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
#include "geometry/create_voronoi_mesh.h"

#if 0
void test_create_single_tetrahedron(void** state)
{
  // Set up points that will yield a single tetrahedral Voronoi cell
  // with vertices (1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1).
  const real_t sqrt2 = sqrt(2.0);
  const real_t sqrt3 = sqrt(3.0);
  const real_t sqrt6 = sqrt(6.0);
  const real_t cos60 = cos(M_PI/3.0);
  const real_t cos30 = cos(M_PI/6.0);
  const real_t cos120 = cos(2.0*M_PI/3.0);

  // Geometric properties.
  real_t side = sqrt2;
  real_t height = sqrt2 * side / sqrt3;
  real_t base_area = 0.25 * sqrt3 * side * side;
  real_t volume = sqrt2 * pow(side, 3) / 12.0;

  // Vectors connecting the center to the various faces.
  real_t d = 0.5 * side / sqrt6;
  vector_t v1 = {.x = 0.0, .y = 0.0,        .z = -d};        // bottom face
  vector_t v2 = {.x = 0.0, .y = -d * cos30, .z = d * cos60}; // -y face
  vector_t v3 = {.x = -d * cos30, .y = d * cos120, .z = d * cos60}; // "-x" face
  vector_t v4 = {.x =  d * cos30, .y = d * cos120, .z = d * cos60}; // "+x" face
  point_t generators[5] = {{.x = 0.0,  .y =  0.0, .z =  0.0},
                           {.x = 2.0*v1.x, .y =  2.0*v1.y, .z = 2.0*v1.z},
                           {.x = 2.0*v2.x, .y =  2.0*v2.y, .z = 2.0*v2.z},
                           {.x = 2.0*v3.x, .y =  2.0*v3.y, .z = 2.0*v3.z},
                           {.x = 2.0*v4.x, .y =  2.0*v4.y, .z = 2.0*v4.z}};

  // Now generate the mesh.
  mesh_t* mesh = create_voronoi_mesh(generators, 5, NULL, 0);
  mesh_verify(mesh);
  assert_int_equal(1, mesh->num_cells);
  assert_int_equal(4, mesh->num_faces);
  assert_int_equal(6, mesh->num_edges);
  assert_int_equal(4, mesh->num_nodes);

  // Verify the volume of the cell.
printf("V = %g vs %g\n", mesh->cells[0].volume, volume);
  assert_true(fabs(volume - mesh->cells[0].volume) < 1e-12);

  // Verify its center position.
  assert_true(fabs(mesh->cells[0].center.x) < 1e-12);
  assert_true(fabs(mesh->cells[0].center.y) < 1e-12);
  assert_true(fabs(mesh->cells[0].center.z) < 1e-12);

  // Verify the areas and centers of its faces.
  for (int f = 0; f < 4; ++f)
  {
    assert_true(fabs(base_area - mesh->faces[f].area) < 1e-12);
  }

  mesh_free(mesh);
}
#endif

void test_create_single_cube(void** state)
{
  // Set up points that will yield a single cubic Voronoi cell.
  point_t generators[7] = {{.x =  0.0, .y =  0.0, .z =  0.0},
                           {.x = -1.0, .y =  0.0, .z =  0.0},
                           {.x =  1.0, .y =  0.0, .z =  0.0},
                           {.x =  0.0, .y = -1.0, .z =  0.0},
                           {.x =  0.0, .y =  1.0, .z =  0.0},
                           {.x =  0.0, .y =  0.0, .z = -1.0},
                           {.x =  0.0, .y =  0.0, .z =  1.0}};

  // Now generate the mesh.
  int_slist_t* deleted_generators = int_slist_new();
  mesh_t* mesh = create_voronoi_mesh(generators, 7, NULL, 0, deleted_generators);
  mesh_verify(mesh);
  assert_int_equal(6, deleted_generators->size);
  assert_int_equal(1, mesh->num_cells);
  assert_int_equal(6, mesh->num_faces);
  assert_int_equal(12, mesh->num_edges);
  assert_int_equal(8, mesh->num_nodes);

  // Verify the volume of the cell.
  assert_true(fabs(1.0 - mesh->cells[0].volume) < 1e-12);

  // Verify its center position.
  assert_true(fabs(mesh->cells[0].center.x) < 1e-12);
  assert_true(fabs(mesh->cells[0].center.y) < 1e-12);
  assert_true(fabs(mesh->cells[0].center.z) < 1e-12);

  // Verify the areas and centers of its faces.
  for (int f = 0; f < 4; ++f)
  {
    real_t face_area = vector_mag(&mesh->faces[f].normal);
    assert_true(fabs(1.0 - face_area) < 1e-12);
  }

  int_slist_free(deleted_generators);
  mesh_free(mesh);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  set_log_level(LOG_DEBUG);
  const UnitTest tests[] = 
  {
    unit_test(test_create_single_cube)
  //  unit_test(test_create_single_tetrahedron)
  };
  return run_tests(tests);
}
