#include "core/unordered_set.h"
#include "io/generate_cell_node_conn.h"

#ifdef __cplusplus
extern "C" {
#endif

void generate_cell_node_conn(mesh_t* mesh,
                             int* face_nodes,
                             int* face_node_offsets,
                             int** cell_nodes,
                             int* cell_node_offsets)
{
  // Make a list of all the cell nodes.
  int_unordered_set_t* all_cell_nodes = int_unordered_set_new();
  for (int c = 0; c < mesh->num_cells; ++c)
  {
    for (int f = 0; f < mesh->cells[c].num_faces; ++f)
    {
      for (int e = 0; e < mesh->cells[c].faces[f]->num_edges; ++e)
      {
        int_unordered_set_insert(all_cell_nodes, mesh->cells[c].faces[f]->edges[e]->node1);
        int_unordered_set_insert(all_cell_nodes, mesh->cells[c].faces[f]->edges[e]->node2);
      }
    }
  }
  *cell_nodes = malloc(all_cell_nodes->size*sizeof(int));
  int_unordered_set_copy_out(all_cell_nodes, *cell_nodes);
}

#ifdef __cplusplus
}
#endif

