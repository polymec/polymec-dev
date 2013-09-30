// Copyright 2012-2013 Jeffrey Johnson.
// 
// This file is part of Polymec, and is licensed under the Apache License, 
// Version 2.0 (the "License"); you may not use this file except in 
// compliance with the License. You may may find the text of the license in 
// the LICENSE file at the top-level source directory, or obtain a copy of 
// it at
// 
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

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
  const double sqrt2 = sqrt(2.0);
  const double sqrt3 = sqrt(3.0);
  const double sqrt6 = sqrt(6.0);
  const double cos60 = cos(M_PI/3.0);
  const double cos30 = cos(M_PI/6.0);
  const double cos120 = cos(2.0*M_PI/3.0);

  // Geometric properties.
  double side = sqrt2;
  double height = sqrt2 * side / sqrt3;
  double base_area = 0.25 * sqrt3 * side * side;
  double volume = sqrt2 * pow(side, 3) / 12.0;

  // Vectors connecting the center to the various faces.
  double d = 0.5 * side / sqrt6;
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
    double face_area = vector_mag(&mesh->faces[f].normal);
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
