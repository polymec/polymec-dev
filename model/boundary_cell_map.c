// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "model/boundary_cell_map.h"

// Constructor for the boundary cell.
static boundary_cell_t* create_boundary_cell(mesh_t* mesh, int cell)
{
  boundary_cell_t* bcell = polymec_malloc(sizeof(boundary_cell_t));
  int num_cell_faces = mesh_cell_num_faces(mesh, cell); 
  bcell->neighbor_cells = polymec_malloc(sizeof(int)*num_cell_faces);
  bcell->boundary_faces = polymec_malloc(sizeof(int)*num_cell_faces);
  bcell->bc_for_face = polymec_malloc(sizeof(void*)*num_cell_faces);
  bcell->opp_faces = polymec_malloc(sizeof(int)*num_cell_faces);
  int pos = 0, face;
  while (mesh_cell_next_face(mesh, cell, &pos, &face))
  {
    int f = pos - 1;
    bcell->neighbor_cells[f] = mesh_face_opp_cell(mesh, face, cell);
    bcell->boundary_faces[f] = (bcell->neighbor_cells[f] == -1) ? face : -1;
    bcell->bc_for_face[f] = NULL;
    bcell->opp_faces[f] = -1;
  }
  return bcell;
}

static void destroy_boundary_cell_entry(int key, boundary_cell_t* cell)
{
  polymec_free(cell->neighbor_cells);
  polymec_free(cell->boundary_faces);
  polymec_free(cell->bc_for_face);
  polymec_free(cell->opp_faces);
  polymec_free(cell);
}

static void destroy_pmap_entry(char* key, void* periodic_map)
{
  int_int_unordered_map_t* pmap = (void*)periodic_map;
  int_int_unordered_map_free(pmap);
}

static string_ptr_unordered_map_t* generate_periodic_maps(mesh_t* mesh, string_ptr_unordered_map_t* bcs)
{
  string_ptr_unordered_map_t* periodic_maps = string_ptr_unordered_map_new();

  // Traverse the boundary conditions and look for periodic_bc objects.
  int pos = 0;
  char* tag;
  void* bc;
  while (string_ptr_unordered_map_next(bcs, &pos, &tag, &bc))
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
      pmap = (int_int_unordered_map_t**)string_ptr_unordered_map_get(periodic_maps, tag2);
    else
      pmap = (int_int_unordered_map_t**)string_ptr_unordered_map_get(periodic_maps, tag1);
    if (pmap != NULL)
    {
      // No destructor for this map, since it's already handled elsewhere.
      string_ptr_unordered_map_insert(periodic_maps, tag, *pmap);
      continue;
    }

    // If we're here, we need to generate the periodic mapping.
    int_int_unordered_map_t* new_pmap = periodic_bc_generate_map(pbc, mesh);

    // Insert it into our table with this tag as the key.
    string_ptr_unordered_map_insert_with_kv_dtor(periodic_maps, tag, new_pmap, destroy_pmap_entry);
  }

  return periodic_maps;
}

boundary_cell_map_t* boundary_cell_map_from_mesh_and_bcs(mesh_t* mesh, string_ptr_unordered_map_t* bcs)
{
  ASSERT(mesh != NULL);
  ASSERT(bcs != NULL);

  // Firstly, we search through our table of boundary conditions and 
  // generate any face-face mappings needed by periodic boundary conditions.
  string_ptr_unordered_map_t* periodic_maps = generate_periodic_maps(mesh, bcs);

  boundary_cell_map_t* boundary_cells = boundary_cell_map_new();

  int pos = 0;
  char* tag;
  void* bc;
  while (string_ptr_unordered_map_next(bcs, &pos, &tag, &bc))
  {
    // Retrieve the tag for this boundary condition.
    ASSERT(mesh_has_tag(mesh->face_tags, tag));
    int num_faces;
    int* faces = mesh_tag(mesh->face_tags, tag, &num_faces);

    // Now create an entry for each boundary cell.
    for (int f = 0; f < num_faces; ++f)
    {
      int face = faces[f];
      if (mesh->face_cells[2*face+1] != -1)
        polymec_error("Boundary face tag '%s' contains an interior face (%d).", tag, face);

      // Get the cell for this boundary face. This is the "boundary" cell, 
      // and must be on the interior of this domain.
      int bcell = mesh->face_cells[2*face];
      ASSERT(bcell < mesh->num_cells);

      boundary_cell_t* boundary_cell;
      boundary_cell_t** bcell_ptr = (boundary_cell_t**)boundary_cell_map_get(boundary_cells, bcell);
      if (bcell_ptr == NULL)
      {
        // Create an entry for this boundary cell.
        boundary_cell = create_boundary_cell(mesh, bcell);
        boundary_cell_map_insert_with_kv_dtor(boundary_cells, bcell, boundary_cell, destroy_boundary_cell_entry);
      }
      else
      {
        // Otherwise, retrieve the one we've already made.
        boundary_cell = *bcell_ptr;
      }

      // Find the local index for this face within the boundary cell.
      int local_face_index = 0;
      int num_cell_faces = mesh->cell_face_offsets[bcell+1] - mesh->cell_face_offsets[bcell];
      for (local_face_index = 0; local_face_index < num_cell_faces; ++local_face_index)
      {
        if (boundary_cell->boundary_faces[local_face_index] == face)
          break;
      }
      ASSERT(local_face_index < num_cell_faces);

      // If the boundary condition is periodic, we have to identify this 
      // boundary face with its other face, and fill in the appropriate 
      // neighbor_cells and opp_faces entries.
      if (pointer_is_periodic_bc(bc))
      {
        int_int_unordered_map_t* periodic_map = *string_ptr_unordered_map_get(periodic_maps, tag);
        int other_face = *int_int_unordered_map_get(periodic_map, face);
        int other_cell = mesh->face_cells[2*other_face];
        ASSERT(other_cell != -1);
        boundary_cell->neighbor_cells[local_face_index] = other_cell;
        boundary_cell->opp_faces[local_face_index] = other_face;
      }
      else
      {
        // Associate the BC.
        boundary_cell->bc_for_face[local_face_index] = bc;
      }
    }
  }

  // Clean up the periodic mappings.
  string_ptr_unordered_map_free(periodic_maps);

  return boundary_cells;
}

