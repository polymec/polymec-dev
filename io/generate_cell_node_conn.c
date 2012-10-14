#include "core/slist.h"
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
  int_slist_t* all_cell_nodes_list = int_slist_new();

}

#ifdef __cplusplus
}
#endif

