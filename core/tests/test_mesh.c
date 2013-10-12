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
#include "core/mesh.h"

void test_single_cell_mesh_no_topo(void** state)
{
  // Create a single hexahedron without topology.
  mesh_t* mesh = mesh_new(1, 0, 6, 12, 8);
  assert_int_equal(1, mesh->num_cells);
  assert_int_equal(0, mesh->num_ghost_cells);
  assert_int_equal(6, mesh->num_faces);
  assert_int_equal(12, mesh->num_edges);
  assert_int_equal(8, mesh->num_nodes);
  mesh_free(mesh);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_single_cell_mesh_no_topo),
  };
  return run_tests(tests);
}
