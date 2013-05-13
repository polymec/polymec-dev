#include "core/unordered_set.h"
#include "io/generate_cell_node_conn.h"

void generate_cell_node_conn(mesh_t* mesh,
                             int* face_nodes,
                             int* face_node_offsets,
                             int** cell_nodes,
                             int* cell_node_offsets)
{
  // Make a list of all the cell nodes.
  int_unordered_set_t* all_cell_nodes = int_unordered_set_new();
  cell_node_offsets[0] = 0;
  for (int c = 0; c < mesh->num_cells; ++c)
  {
    for (int f = 0; f < mesh->cells[c].num_faces; ++f)
    {
      int face_id = mesh->cells[c].faces[f] - &mesh->faces[0];
      ASSERT(face_id < mesh->num_faces);
      int offset = face_node_offsets[face_id];
      int nn = face_node_offsets[face_id+1] - offset;
      for (int n = 0; n < nn; ++n)
        int_unordered_set_insert(all_cell_nodes, face_nodes[offset+n]);
    }
    cell_node_offsets[c+1] = all_cell_nodes->size;
  }
  *cell_nodes = malloc(all_cell_nodes->size*sizeof(int));
  int pos = 0, node, i = 0;
  while (int_unordered_set_next(all_cell_nodes, &pos, &node))
    (*cell_nodes)[i++] = node;
  int_unordered_set_free(all_cell_nodes);
}

