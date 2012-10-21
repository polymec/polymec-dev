#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmockery.h"
#include "cubic_lattice.h"
#include "create_cubic_lattice_mesh.h"
#include "vtk_plot_io.h"

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
  // Create a 10x10x10 cubic lattice mesh.
  mesh_t* mesh = create_cubic_lattice_mesh(10, 10, 10);

  // Plot it.
  io_interface_t* plot = vtk_plot_io_new(MPI_COMM_SELF, 0, false);
  io_open(plot, "cubic_lattice_10x10x10", ".", IO_WRITE);
  io_dataset_t* dataset = io_dataset_new("default", 1, 0);
  io_dataset_write_mesh(dataset, mesh);
  double ones[1000];
  for (int c = 0; c < 1000; ++c)
    ones[c] = 1.0;
  io_dataset_write_field(dataset, "solution", ones, 1000, MESH_CELL);
  io_append_dataset(plot, dataset);
  io_close(plot);

  // Clean up.
  io_free(plot);
  mesh_free(mesh);
}

int main(int argc, char* argv[]) 
{
  arbi_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_create_cubic_lattice_mesh),
    unit_test(test_plot_cubic_lattice_mesh)
  };
  return run_tests(tests);
}
