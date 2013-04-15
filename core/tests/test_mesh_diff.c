#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmockery.h"
#include "core/mesh.h"
#include "core/mesh_diff.h"

void test_single_cell_mesh_no_topo(void** state)
{
  mesh_t* mesh = mesh_new(0, 0, 0, 0, 0);

  // Create a single hexahedron without topology.
  mesh_diff_t* diff = mesh_diff_new();

  // Nodes.
  point_t nodes[8] = {{.x = 0.0, .y = 0.0, .z = 0.0}, 
                      {.x = 1.0, .y = 0.0, .z = 0.0},
                      {.x = 1.0, .y = 1.0, .z = 0.0},
                      {.x = 0.0, .y = 1.0, .z = 0.0},
                      {.x = 0.0, .y = 0.0, .z = 1.0}, 
                      {.x = 1.0, .y = 0.0, .z = 1.0},
                      {.x = 1.0, .y = 1.0, .z = 1.0},
                      {.x = 0.0, .y = 1.0, .z = 1.0}};
  for (int n = 0; n < 8; ++n)
  {
    mesh_delta_t* append_node = append_node_mesh_delta_new(&nodes[n]);
    mesh_diff_append(diff, append_node);
  }

  // Edges.
  for (int e = 0; e < 12; ++e)
  {
    mesh_delta_t* append_edge = append_mesh_delta_new(MESH_EDGE);
    mesh_diff_append(diff, append_edge);
  }

  // Faces.
  for (int f = 0; f < 6; ++f)
  {
    mesh_delta_t* append_face = append_mesh_delta_new(MESH_FACE);
    mesh_diff_append(diff, append_face);
  }

  // Cell.
  mesh_delta_t* append_cell = append_mesh_delta_new(MESH_CELL);
  mesh_diff_append(diff, append_cell);

  // Apply the diff.
  mesh_diff_apply(diff, mesh);

  // Check the mesh.
  assert_int_equal(1, mesh->num_cells);
  assert_int_equal(0, mesh->num_ghost_cells);
  assert_int_equal(6, mesh->num_faces);
  assert_int_equal(12, mesh->num_edges);
  assert_int_equal(8, mesh->num_nodes);

  // Roll back the diff.
  assert_true(mesh_property(mesh, "last_diff") != NULL);
  mesh_diff_rollback(diff, mesh);

  // Check the mesh.
  assert_int_equal(0, mesh->num_cells);
  assert_int_equal(0, mesh->num_ghost_cells);
  assert_int_equal(0, mesh->num_faces);
  assert_int_equal(0, mesh->num_edges);
  assert_int_equal(0, mesh->num_nodes);

  // Clean up.
  mesh_diff_free(diff);
  mesh_free(mesh);
}

void test_single_cell_mesh(void** state)
{
  mesh_t* mesh = mesh_new(0, 0, 0, 0, 0);

  // Create a single hexahedron.
  mesh_diff_t* diff = mesh_diff_new();

  // Nodes.
  point_t nodes[8] = {{.x = 0.0, .y = 0.0, .z = 0.0}, 
                      {.x = 1.0, .y = 0.0, .z = 0.0},
                      {.x = 1.0, .y = 1.0, .z = 0.0},
                      {.x = 0.0, .y = 1.0, .z = 0.0},
                      {.x = 0.0, .y = 0.0, .z = 1.0}, 
                      {.x = 1.0, .y = 0.0, .z = 1.0},
                      {.x = 1.0, .y = 1.0, .z = 1.0},
                      {.x = 0.0, .y = 1.0, .z = 1.0}};
  for (int n = 0; n < 8; ++n)
  {
    mesh_delta_t* append_node = append_node_mesh_delta_new(&nodes[n]);
    mesh_diff_append(diff, append_node);
  }

  // Edges.
  int edges[12][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}, 
                      {4, 5}, {5, 6}, {6, 7}, {7, 4},
                      {0, 4}, {1, 5}, {2, 6}, {3, 7}};
  for (int e = 0; e < 12; ++e)
  {
    mesh_delta_t* append_edge = append_mesh_delta_new(MESH_EDGE);
    mesh_diff_append(diff, append_edge);

    mesh_delta_t* attach_1 = attach_mesh_delta_new(MESH_NODE, edges[e][0], e);
    mesh_diff_append(diff, attach_1);
    mesh_delta_t* attach_2 = attach_mesh_delta_new(MESH_NODE, edges[e][1], e);
    mesh_diff_append(diff, attach_2);
  }

  // Faces.
  int faces[6][4] = {{0, 9, 5, 8}, {1, 10, 6, 9}, 
                     {2, 11, 7, 10}, {3, 8, 7, 11},
                     {0, 1, 2, 3}, {4, 5, 6, 7}};
  for (int f = 0; f < 6; ++f)
  {
    mesh_delta_t* append_face = append_mesh_delta_new(MESH_FACE);
    mesh_diff_append(diff, append_face);

    mesh_delta_t* attach_1 = attach_mesh_delta_new(MESH_EDGE, faces[f][0], f);
    mesh_diff_append(diff, attach_1);
    mesh_delta_t* attach_2 = attach_mesh_delta_new(MESH_EDGE, faces[f][1], f);
    mesh_diff_append(diff, attach_2);
    mesh_delta_t* attach_3 = attach_mesh_delta_new(MESH_EDGE, faces[f][2], f);
    mesh_diff_append(diff, attach_3);
    mesh_delta_t* attach_4 = attach_mesh_delta_new(MESH_EDGE, faces[f][3], f);
    mesh_diff_append(diff, attach_4);
  }

  // Cell.
  mesh_delta_t* append_cell = append_mesh_delta_new(MESH_CELL);
  mesh_diff_append(diff, append_cell);

  mesh_delta_t* attach_1 = attach_mesh_delta_new(MESH_FACE, 0, 0);
  mesh_diff_append(diff, attach_1);
  mesh_delta_t* attach_2 = attach_mesh_delta_new(MESH_FACE, 1, 0);
  mesh_diff_append(diff, attach_2);
  mesh_delta_t* attach_3 = attach_mesh_delta_new(MESH_FACE, 2, 0);
  mesh_diff_append(diff, attach_3);
  mesh_delta_t* attach_4 = attach_mesh_delta_new(MESH_FACE, 3, 0);
  mesh_diff_append(diff, attach_4);
  mesh_delta_t* attach_5 = attach_mesh_delta_new(MESH_FACE, 4, 0);
  mesh_diff_append(diff, attach_5);
  mesh_delta_t* attach_6 = attach_mesh_delta_new(MESH_FACE, 5, 0);
  mesh_diff_append(diff, attach_6);

  // Apply the diff.
  mesh_diff_apply(diff, mesh);

  // Check the mesh.
  assert_int_equal(1, mesh->num_cells);
  assert_int_equal(6, mesh->cells[0].num_faces);
  assert_int_equal(0, mesh->num_ghost_cells);
  assert_int_equal(6, mesh->num_faces);
  assert_int_equal(12, mesh->num_edges);
  assert_int_equal(8, mesh->num_nodes);
  for (int f = 0; f < 6; ++f)
  {
    assert_int_equal(4, mesh->cells[0].faces[f]->num_edges);
    for (int e = 0; e < 4; ++e)
    {
      assert_int_equal(faces[f][e], mesh->cells[0].faces[f]->edges[e] - mesh->edges);
      assert_int_equal(edges[e][0], mesh->cells[0].faces[f]->edges[e]->node1 - mesh->nodes);
      assert_int_equal(edges[e][1], mesh->cells[0].faces[f]->edges[e]->node2 - mesh->nodes);
    }
  }

  // Roll back the diff.
  mesh_diff_rollback(diff, mesh);

  // Check the mesh.
  assert_int_equal(0, mesh->num_cells);
  assert_int_equal(0, mesh->num_ghost_cells);
  assert_int_equal(0, mesh->num_faces);
  assert_int_equal(0, mesh->num_edges);
  assert_int_equal(0, mesh->num_nodes);

  // Clean up.
  mesh_diff_free(diff);
  mesh_free(mesh);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_single_cell_mesh_no_topo),
    unit_test(test_single_cell_mesh)
  };
  return run_tests(tests);
}
