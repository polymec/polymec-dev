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
#include "core/write_silo.h"
#include "geometry/create_cubic_lattice_mesh.h"

#undef HAVE_HDF5
#undef HAVE_TETGEN
#include "polytope_c.h"

void test_create_cubic_lattice_mesh(void** state)
{
  // Create a 10x10x10 cubic lattice mesh.
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  mesh_t* mesh = create_cubic_lattice_mesh(10, 10, 10, &bbox);
  mesh_verify(mesh);
  assert_int_equal(10*10*10, mesh->num_cells);
  assert_int_equal(0, mesh->num_ghost_cells);
  assert_int_equal(mesh->num_faces, 10*11*11);
  assert_int_equal(mesh->num_edges, 10*10*11);
  assert_int_equal(11*11*11, mesh->num_nodes);
  mesh_free(mesh);
}

void test_plot_cubic_lattice_mesh(void** state)
{
  // Create a 4x4x4 cubic lattice mesh.
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  mesh_t* mesh = create_cubic_lattice_mesh(4, 4, 4, &bbox);

  // Plot it.
  double ones[4*4*4];
  for (int c = 0; c < 4*4*4; ++c)
    ones[c] = 1.0*c;
  string_ptr_unordered_map_t* fields = string_ptr_unordered_map_new();
  string_ptr_unordered_map_insert(fields, "solution", ones);
  write_silo(mesh, NULL, NULL, NULL, fields, "cubic_lattice_4x4x4", ".",
             0, 0.0, MPI_COMM_SELF, 1, 0);

  // Clean up.
  string_ptr_unordered_map_free(fields);
  mesh_free(mesh);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_create_cubic_lattice_mesh),
    unit_test(test_plot_cubic_lattice_mesh)
  };
  return run_tests(tests);
}
