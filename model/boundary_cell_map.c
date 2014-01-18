// Copyright (c) 2012-2014, Jeffrey N. Johnson
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this 
// list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice, 
// this list of conditions and the following disclaimer in the documentation 
// and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "model/boundary_cell_map.h"

// Constructor for the boundary cell.
static boundary_cell_t* create_boundary_cell()
{
  boundary_cell_t* bcell = malloc(sizeof(boundary_cell_t));
  bcell->neighbor_cells = NULL;
  bcell->num_neighbor_cells = 0;
  bcell->boundary_faces = NULL;
  bcell->num_boundary_faces = 0;
  bcell->bc_for_face = NULL;
  bcell->opp_faces = NULL;
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
  if (cell->opp_faces != NULL)
    free(cell->opp_faces);
  free(cell);
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

    // Now create an entry for each boundary cell and count boundary
    // faces and neighbors.
    for (int f = 0; f < num_faces; ++f)
    {
      int face = faces[f];
      ASSERT(mesh->face_cells[2*face+1] == -1); // ... true for now.

      // Get the cell for this boundary face. This is the boundary cell.
      int bcell = mesh->face_cells[2*face];

      // Have we created an entry for this one yet? If not, do so.
      boundary_cell_t* boundary_cell;
      if (!boundary_cell_map_contains(boundary_cells, bcell))
      {
        boundary_cell = create_boundary_cell();
        boundary_cell_map_insert_with_kv_dtor(boundary_cells, bcell, boundary_cell, destroy_boundary_cell_entry);

        // Gather the interior faces for the cell.
        int pos = 0, ff;
        while (mesh_cell_next_face(mesh, bcell, &pos, &ff))
        {
          if (mesh_face_opp_cell(mesh, bcell, ff) != -1)
            boundary_cell->num_neighbor_cells++;
        }
        boundary_cell->neighbor_cells = malloc(sizeof(int)*boundary_cell->num_neighbor_cells);
        int nc = 0;
        pos = 0;
        while (mesh_cell_next_face(mesh, bcell, &pos, &ff))
        {
          int opp_cell = mesh_face_opp_cell(mesh, ff, bcell);
          if (opp_cell != -1)
            boundary_cell->neighbor_cells[nc++] = opp_cell;
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
      bcell->opp_faces = malloc(sizeof(int)*bcell->num_boundary_faces);
    }
  }

  // Now go back through and set the boundary faces and boundary 
  // conditions for each cell.
  pos = 0;
  while (string_ptr_unordered_map_next(bcs, &pos, &tag, &bc))
  {
    // Retrieve the tag for this boundary condition.
    ASSERT(mesh_has_tag(mesh->face_tags, tag));
    int num_faces;
    int* faces = mesh_tag(mesh->face_tags, tag, &num_faces);

    // Now create an entry for each boundary cell and count boundary
    // faces and neighbors.
    for (int f = 0; f < num_faces; ++f)
    {
      int face = faces[f];
      int bcell = mesh->face_cells[2*face];

      boundary_cell_t* boundary_cell = *boundary_cell_map_get(boundary_cells, bcell);
      ASSERT(boundary_cell != NULL);

      int i = 0;
      while (boundary_cell->boundary_faces[i] != -1) ++i;
      boundary_cell->boundary_faces[i] = face;

      // If the boundary condition is periodic, we have to identify this 
      // boundary face with its other face, and fill in the appropriate 
      // opp_cell entry.
      if (pointer_is_periodic_bc(bc))
      {
        boundary_cell->bc_for_face[i] = NULL;
        int_int_unordered_map_t* periodic_map = *string_ptr_unordered_map_get(periodic_maps, tag);
        int other_face_index = *int_int_unordered_map_get(periodic_map, face);
        boundary_cell->opp_faces[i] = other_face_index;
      }
      else
      {
        boundary_cell->bc_for_face[i] = bc;
        boundary_cell->opp_faces[i] = -1;
      }
    }
  }

  // Clean up the periodic mappings.
  string_ptr_unordered_map_free(periodic_maps);

  return boundary_cells;
}

