#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmockery.h"
#include "mesh.h"
#include "edit_mesh.h"

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

void test_single_cell_mesh(void** state)
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
                          {5, 1, 6, 10},
                          {6, 2, 7, 11},
                          {8, 9, 10, 11}};
  int edge_nodes[12][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0},
                           {0, 4}, {1, 5}, {2, 6}, {3, 7},
                           {4, 5}, {5, 6}, {6, 7}, {7, 0}};

  for (int f = 0; f < 6; ++f)
  {
    mesh_add_face_to_cell(mesh, &mesh->faces[f], &mesh->cells[0]);
    for (int e = 0; e < 4; ++e)
    {
      mesh_add_edge_to_face(mesh, &mesh->edges[face_edges[f][e]], &mesh->faces[f]);
      int fe = face_edges[f][e];
      for (int n = 0; n < 2; ++n)
        mesh->edges[fe].node1 = &mesh->nodes[edge_nodes[fe][n]];
    }
  }
  mesh_free(mesh);
}

int main(int argc, char* argv[]) 
{
  arbi_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_single_cell_mesh_no_topo),
    unit_test(test_single_cell_mesh)
  };
  return run_tests(tests);
}
