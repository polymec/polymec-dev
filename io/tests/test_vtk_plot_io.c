#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmockery.h"
#include "core/edit_mesh.h"
#include "io/vtk_plot_io.h"

void test_plot_single_cell_mesh(void** state)
{
  // Create a single hexahedron.
  mesh_t* mesh = mesh_new(1, 0, 6, 12, 8);
  assert_int_equal(1, mesh->num_cells);
  assert_int_equal(0, mesh->num_ghost_cells);
  assert_int_equal(6, mesh->num_faces);
  assert_int_equal(12, mesh->num_edges);
  assert_int_equal(8, mesh->num_nodes);

  int face_edges[6][4] = {{0, 1, 2, 3}, 
                          {3, 4, 11, 7},
                          {4, 0, 5, 8},
                          {5, 1, 6, 9},
                          {6, 2, 7, 10},
                          {8, 9, 10, 11}};
  int edge_nodes[12][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0},
                           {0, 4}, {1, 5}, {2, 6}, {3, 7},
                           {4, 5}, {5, 6}, {6, 7}, {7, 4}};
  node_t nodes[8] = {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0},
                     {1.0, 1.0, 0.0}, {0.0, 1.0, 0.0},
                     {0.0, 0.0, 1.0}, {1.0, 0.0, 1.0},
                     {1.0, 1.0, 1.0}, {0.0, 1.0, 1.0}};

  for (int f = 0; f < 6; ++f)
  {
    mesh_attach_face_to_cell(mesh, &mesh->faces[f], &mesh->cells[0]);
    for (int e = 0; e < 4; ++e)
    {
      mesh_attach_edge_to_face(mesh, &mesh->edges[face_edges[f][e]], &mesh->faces[f]);
      int fe = face_edges[f][e];
      mesh->edges[fe].node1 = &mesh->nodes[edge_nodes[fe][0]];
      mesh->edges[fe].node2 = &mesh->nodes[edge_nodes[fe][1]];
    }
  }

  assert_int_equal(6, mesh->cells[0].num_faces);
  for (int f = 0; f < 6; ++f)
    assert_int_equal(4, mesh->faces[f].num_edges);

  for (int n = 0; n < 8; ++n)
    mesh->nodes[n] = nodes[n];

  // Plot it.
  io_interface_t* plot = vtk_plot_io_new(MPI_COMM_SELF, 0, false);
  io_open(plot, "hex", ".", IO_WRITE);
  io_dataset_t* dataset = io_dataset_new("default");
  io_dataset_put_mesh(dataset, mesh);
  double one = 1.0;
  io_dataset_put_field(dataset, "solution", &one, 1, MESH_CELL, true);
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
    unit_test(test_plot_single_cell_mesh)
  };
  return run_tests(tests);
}
