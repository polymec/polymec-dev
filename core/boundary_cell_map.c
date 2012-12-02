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
  bcell->opp_cells = NULL;
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
  if (cell->opp_cells != NULL)
    free(cell->opp_cells);
  free(cell);
}

static void destroy_pmap_entry(char* key, void* periodic_map)
{
  int_int_unordered_map_t* pmap = (void*)periodic_map;
  int_int_unordered_map_free(pmap);
}

static str_ptr_unordered_map_t* generate_periodic_maps(mesh_t* mesh, str_ptr_unordered_map_t* bcs)
{
  str_ptr_unordered_map_t* periodic_maps = str_ptr_unordered_map_new();

  // Traverse the boundary conditions and look for periodic_bc objects.
  int pos = 0;
  char* tag;
  void* bc;
  while (str_ptr_unordered_map_next(bcs, &pos, &tag, &bc))
  {
    periodic_bc_t* pbc = (periodic_bc_t*)bc;
    if (!periodic_bc_is_valid(pbc)) continue;

    // Make sure that the periodic BC has tags that are valid.
    char *tag1, *tag2;
    periodic_bc_get_tags(pbc, &tag1, &tag2);
    if (!mesh_has_tag(mesh->face_tags, (const char*)tag1))
      polymec_error("Face tag '%s' for periodic BC not found in the mesh.", tag1);
    if (!mesh_has_tag(mesh->face_tags, (const char*)tag2))
      polymec_error("Face tag '%s' for periodic BC not found in the mesh.", tag2);

    // If a periodic map has already been generated for this boundary
    // condition, assign it and move on.
    int_int_unordered_map_t** pmap;
    if (!strcmp(tag, tag1))
      pmap = (int_int_unordered_map_t**)str_ptr_unordered_map_get(periodic_maps, tag2);
    else
      pmap = (int_int_unordered_map_t**)str_ptr_unordered_map_get(periodic_maps, tag1);
    if (pmap != NULL)
    {
      // No destructor for this map, since it's already handled elsewhere.
      str_ptr_unordered_map_insert(periodic_maps, tag, *pmap);
      continue;
    }

    // If we're here, we need to generate the periodic mapping.
    int_int_unordered_map_t* new_pmap = periodic_bc_generate_map(pbc, mesh);

    // Insert it into our table with this tag as the key.
    str_ptr_unordered_map_insert_with_dtor(periodic_maps, tag, new_pmap, destroy_pmap_entry);
  }

  return periodic_maps;
}

boundary_cell_map_t* boundary_cell_map_from_mesh_and_bcs(mesh_t* mesh, str_ptr_unordered_map_t* bcs)
{
  ASSERT(mesh != NULL);
  ASSERT(bcs != NULL);

  // Firstly, we search through our table of boundary conditions and 
  // generate any face-face mappings needed by periodic boundary conditions.
  str_ptr_unordered_map_t* periodic_maps = generate_periodic_maps(mesh, bcs);

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
      bcell->opp_cells = malloc(sizeof(cell_t*)*bcell->num_boundary_faces);
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

      // If the boundary condition is periodic, we have to identify this 
      // boundary face with its other face, and fill in the appropriate 
      // opp_cell entry.
      if (pointer_is_periodic_bc(bc))
      {
        boundary_cell->bc_for_face[i] = NULL;
        int_int_unordered_map_t* periodic_map = *str_ptr_unordered_map_get(periodic_maps, tag);
        int face_index = face - &mesh->faces[0];
        int other_face_index = *int_int_unordered_map_get(periodic_map, face_index);
        face_t* other_face = &mesh->faces[other_face_index];
        boundary_cell->opp_cells[i] = other_face->cell1;
      }
      else
      {
        boundary_cell->bc_for_face[i] = bc;
        boundary_cell->opp_cells[i] = NULL;
      }
    }
  }

  // Clean up the periodic mappings.
  str_ptr_unordered_map_free(periodic_maps);

  return boundary_cells;
}

#ifdef __cplusplus
}
#endif

