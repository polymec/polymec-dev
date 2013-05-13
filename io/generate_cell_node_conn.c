#include "core/slist.h"
#include "io/generate_cell_node_conn.h"

void generate_cell_node_conn(mesh_t* mesh,
                             int* face_nodes,
                             int* face_node_offsets,
                             int** cell_nodes,
                             int* cell_node_offsets)
{
  // Make a list of all the cell nodes.
  int_slist_t* all_cell_nodes = int_slist_new();
  cell_node_offsets[0] = 0;
  for (int c = 0; c < mesh->num_cells; ++c)
  {
    if (mesh->cells[c].num_faces < 4)
      polymec_error("generate_cell_node_conn: Cell %d has fewer than 4 faces.", c);

    for (int f = 0; f < mesh->cells[c].num_faces; ++f)
    {
      int face_id = mesh->cells[c].faces[f] - &mesh->faces[0];
      ASSERT(face_id < mesh->num_faces);
      int offset = face_node_offsets[face_id];
      int nn = face_node_offsets[face_id+1] - offset;
      ASSERT(nn > 0);
      for (int n = 0; n < nn; ++n)
        int_slist_append(all_cell_nodes, face_nodes[offset+n]);
    }
    ASSERT(all_cell_nodes->size > cell_node_offsets[c]);
    cell_node_offsets[c+1] = all_cell_nodes->size;
  }
  *cell_nodes = malloc(all_cell_nodes->size*sizeof(int));
  int i = 0;
  while (!int_slist_empty(all_cell_nodes))
    (*cell_nodes)[i++] = int_slist_pop(all_cell_nodes, NULL);
  int_slist_free(all_cell_nodes);
}

