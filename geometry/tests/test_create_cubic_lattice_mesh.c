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
#include "geometry/cubic_lattice.h"
#include "geometry/create_cubic_lattice_mesh.h"
#include "io/vtk_plot_io.h"

void test_create_cubic_lattice_mesh(void** state)
{
  // Create a 10x10x10 cubic lattice mesh.
  cubic_lattice_t* lattice = cubic_lattice_new(10, 10, 10);
  mesh_t* mesh = create_cubic_lattice_mesh(10, 10, 10);
  mesh_verify(mesh);
  assert_int_equal(mesh->num_cells, cubic_lattice_num_cells(lattice));
  assert_int_equal(1000, mesh->num_cells);
  assert_int_equal(0, mesh->num_ghost_cells);
  assert_int_equal(mesh->num_faces, cubic_lattice_num_faces(lattice));
  assert_int_equal(mesh->num_edges, cubic_lattice_num_edges(lattice));
  assert_int_equal(mesh->num_nodes, cubic_lattice_num_nodes(lattice));
  assert_int_equal(11*11*11, mesh->num_nodes);
  mesh_free(mesh);
}

void test_plot_cubic_lattice_mesh(void** state)
{
  // Create a 4x4x4 cubic lattice mesh.
  mesh_t* mesh = create_cubic_lattice_mesh(4, 4, 4);

  // Plot it.
  io_interface_t* plot = vtk_plot_io_new(MPI_COMM_SELF, 0, false);
  io_open(plot, "cubic_lattice_4x4x4", ".", IO_WRITE);
  io_dataset_t* dataset = io_dataset_new("default");
  io_dataset_put_mesh(dataset, mesh);
  double ones[4*4*4];
  for (int c = 0; c < 4*4*4; ++c)
    ones[c] = 1.0*c;
  io_dataset_put_field(dataset, "solution", ones, 1, MESH_CELL, true);
  io_append_dataset(plot, dataset);
  io_close(plot);

  // Clean up.
  io_free(plot);
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
