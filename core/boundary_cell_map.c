#include "core/boundary_cell_map.h"

#ifdef __cplusplus
extern "C" {
#endif

// Constructor for the boundary cell.
static boundary_cell_t* create_boundary_cell()
{
  boundary_cell_t* bcell = malloc(sizeof(boundary_cell_t));
  bcell->neighbor_cells = NULL;
  bcell->num_neighbor_cells = 0;
  bcell->boundary_faces = NULL;
  bcell->num_boundary_faces = 0;
  bcell->bc_for_face = NULL;
  return bcell;
}

static void destroy_boundary_cell_entry(int key, boundary_cell_t* cell)
{
  if (cell->neighbor_cells != NULL)
    free(cell->neighbor_cells);
  if (cell->boundary_faces != NULL)
    free(cell->boundary_faces);
  if (cell->bc_for_face != NULL)
    free(cell->bc_for_face);
  free(cell);
}

boundary_cell_map_t* boundary_cell_map_from_mesh_and_bcs(mesh_t* mesh, str_ptr_unordered_map_t* bcs)
{
  ASSERT(mesh != NULL);
  ASSERT(bcs != NULL);

  boundary_cell_map_t* boundary_cells = boundary_cell_map_new();

  int pos = 0;
  char* tag;
  void* bc;
  while (str_ptr_unordered_map_next(bcs, &pos, &tag, &bc))
  {
    // Retrieve the tag for this boundary condition.
    ASSERT(mesh_has_tag(mesh->face_tags, tag));
    int num_faces;
    int* faces = mesh_tag(mesh->face_tags, tag, &num_faces);

    // Now create an entry for each boundary cell and count boundary
    // faces and neighbors.
    for (int f = 0; f < num_faces; ++f)
    {
      face_t* face = &mesh->faces[faces[f]];
      ASSERT(face->cell2 == NULL); // ... true for now.

      // Get the cell for this boundary face. This is the boundary cell.
      cell_t* cell = face->cell1;
      int bcell = cell - &mesh->cells[0];

      // Have we created an entry for this one yet? If not, do so.
      boundary_cell_t* boundary_cell;
      if (!boundary_cell_map_contains(boundary_cells, bcell))
      {
        boundary_cell = create_boundary_cell();
        boundary_cell_map_insert_with_dtor(boundary_cells, bcell, boundary_cell, destroy_boundary_cell_entry);

        // Gather the interior faces for the cell.
        for (int ff = 0; ff < cell->num_faces; ++ff)
        {
          if (face_opp_cell(cell->faces[ff], cell) != NULL)
            boundary_cell->num_neighbor_cells++;
        }
        boundary_cell->neighbor_cells = malloc(sizeof(int)*boundary_cell->num_neighbor_cells);
        int nc = 0;
        for (int ff = 0; ff < cell->num_faces; ++ff)
        {
          cell_t* opp_cell = face_opp_cell(cell->faces[ff], cell);
          if (opp_cell != NULL)
          {
            int opp = opp_cell - &mesh->cells[0];
            boundary_cell->neighbor_cells[nc++] = opp;
          }
        }
      }
      else
      {
        // Just retrieve the existing boundary cell.
        boundary_cell = *boundary_cell_map_get(boundary_cells, bcell);
      }

      // Now increment the boundary face count for this cell.
      boundary_cell->num_boundary_faces++;
    }
  }

  // Allocate storage for the boundary faces and BCs within the cells.
  {
    int pos = 0;
    int bcell_index;
    boundary_cell_t* bcell;
    while (boundary_cell_map_next(boundary_cells, &pos, &bcell_index, &bcell))
    {
      bcell->boundary_faces = malloc(sizeof(int)*bcell->num_boundary_faces);
      for (int f = 0; f < bcell->num_boundary_faces; ++f)
        bcell->boundary_faces[f] = -1;
      bcell->bc_for_face = malloc(sizeof(void*)*bcell->num_boundary_faces);
    }
  }

  // Now go back through and set the boundary faces and boundary 
  // conditions for each cell.
  pos = 0;
  while (str_ptr_unordered_map_next(bcs, &pos, &tag, &bc))
  {
    // Retrieve the tag for this boundary condition.
    ASSERT(mesh_has_tag(mesh->face_tags, tag));
    int num_faces;
    int* faces = mesh_tag(mesh->face_tags, tag, &num_faces);

    // Now create an entry for each boundary cell and count boundary
    // faces and neighbors.
    for (int f = 0; f < num_faces; ++f)
    {
      face_t* face = &mesh->faces[faces[f]];
      cell_t* cell = face->cell1;
      int bcell = cell - &mesh->cells[0];

      boundary_cell_t* boundary_cell = *boundary_cell_map_get(boundary_cells, bcell);
      ASSERT(boundary_cell != NULL);

      int i = 0;
      while (boundary_cell->boundary_faces[i] != -1) ++i;
      boundary_cell->boundary_faces[i] = faces[f];
      boundary_cell->bc_for_face[i] = bc;
    }
  }
  return boundary_cells;
}

#ifdef __cplusplus
}
#endif

