#include <stdlib.h>
#include "mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

//------------------------------------------------------------------------
mesh_t* mesh_new(int num_cells, int num_ghost_cells, int num_faces,
                 int num_edges, int num_nodes)
{
  ASSERT(num_cells > 0);
  ASSERT(num_ghost_cells > 0);
  ASSERT(num_faces > 0);
  ASSERT(num_edges > 0);
  ASSERT(num_nodes > 0);
  mesh_t* mesh = malloc(sizeof(mesh_t));

  mesh->cells = malloc(sizeof(cell_t*)*(num_cells+num_ghost_cells));
  for (int i = 0; i < (num_cells + num_ghost_cells); ++i)
    mesh->cells[i] = malloc(sizeof(cell_t));
  mesh->num_cells = num_cells;
  mesh->num_ghost_cells = num_ghost_cells;

  mesh->faces = malloc(sizeof(face_t*)*num_faces);
  for (int i = 0; i < num_faces; ++i)
    mesh->faces[i] = malloc(sizeof(face_t));
  mesh->num_faces = num_faces;

  mesh->edges = malloc(sizeof(edge_t*)*num_edges);
  for (int i = 0; i < num_edges; ++i)
    mesh->edges[i] = malloc(sizeof(edge_t));
  mesh->num_edges = num_edges;

  mesh->nodes = malloc(sizeof(node_t*)*num_nodes);
  for (int i = 0; i < num_nodes; ++i)
    mesh->nodes[i] = malloc(sizeof(node_t));
  mesh->num_nodes = num_nodes;

  return mesh;
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
void mesh_free(mesh_t* mesh)
{
  ASSERT(mesh != NULL);

  for (int i = 0; i < (mesh->num_cells + mesh->num_ghost_cells); ++i)
  {
    if (mesh->cells[i]->faces != NULL)
      free(mesh->cells[i]->faces);
    free(mesh->cells[i]);
  }
  free(mesh->cells);

  for (int i = 0; i < mesh->num_faces; ++i)
  {
    if (mesh->faces[i]->edges != NULL)
      free(mesh->faces[i]->edges);
    free(mesh->faces[i]);
  }
  free(mesh->faces);

  for (int i = 0; i < mesh->num_edges; ++i)
    free(mesh->edges[i]);
  free(mesh->edges);

  for (int i = 0; i < mesh->num_nodes; ++i)
    free(mesh->nodes[i]);
  free(mesh->nodes);

  free(mesh);
}
//------------------------------------------------------------------------

#ifdef __cplusplus
}
#endif

