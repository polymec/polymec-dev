#include <stdlib.h>
#include "lite_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

//------------------------------------------------------------------------
lite_mesh_t* lite_mesh_new(int num_cells, int num_ghost_cells, int num_faces)
{
  ASSERT(num_cells > 0);
  ASSERT(num_ghost_cells > 0);
  ASSERT(num_faces > 0);
  lite_mesh_t* mesh = malloc(sizeof(lite_mesh_t));
  mesh->cells = malloc(sizeof(lite_cell_t*)*(num_cells+num_ghost_cells));
  for (int i = 0; i < (num_cells + num_ghost_cells); ++i)
    mesh->cells[i] = malloc(sizeof(lite_cell_t));
  mesh->num_cells = num_cells;
  mesh->num_ghost_cells = num_ghost_cells;
  mesh->faces = malloc(sizeof(lite_face_t*)*num_faces);
  for (int i = 0; i < num_faces; ++i)
    mesh->faces[i] = malloc(sizeof(lite_face_t));
  mesh->num_faces = num_faces;
  return mesh;
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
void lite_mesh_free(lite_mesh_t* mesh)
{
  ASSERT(mesh != NULL);
  for (int i = 0; i < (mesh->num_cells + mesh->num_ghost_cells); ++i)
    free(mesh->cells[i]);
  free(mesh->cells);
  for (int i = 0; i < mesh->num_faces; ++i)
    free(mesh->faces[i]);
  free(mesh->faces);
  free(mesh);
}
//------------------------------------------------------------------------

#ifdef __cplusplus
}
#endif

