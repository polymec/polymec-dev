// Copyright (c) 2012-2015, Jeffrey N. Johnson
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

#ifndef POLYMEC_BOUNDARY_CELL_MAP_H
#define POLYMEC_BOUNDARY_CELL_MAP_H

#include "core/mesh.h"
#include "core/unordered_map.h"
#include "model/periodic_bc.h"

// This data structure provides metadata for cells that touch the boundary 
// of a domain.
typedef struct
{
  // Cells neighboring this boundary cell. This array is equal in length to /
  // the number of faces for the cell, and neighbor_cells[i] is equal to 
  // the index of the neighboring cell, or -1 if no such cell exists. This 
  // array will be different from the mesh's connectivity data in the case 
  // of periodic boundary conditions.
  int* neighbor_cells;

  // Boundary faces attached to this cell. boundary_faces[i] contains the 
  // index of the ith boundary face within the cell, or -1 if the face is 
  // an interior face.
  int* boundary_faces;

  // bc_for_face[i] stores a pointer to the boundary condition for the 
  // ith face of this cell, or NULL if the ith face is not a boundary face, 
  // OR if the ith face has a periodic boundary condition.
  void** bc_for_face;

  // For faces with periodic boundary conditions, this array has non-negative
  // entries identifying the face on the opposite side of the ith face.
  // opp_faces[i] is -1 if the ith boundary face does not have a periodic 
  // boundary condition. 
  int* opp_faces;
} boundary_cell_t;

// Define the mapping type.
DEFINE_UNORDERED_MAP(boundary_cell_map, int, boundary_cell_t*, int_hash, int_equals)

// Given a mesh and a table of boundary conditions (mesh tags mapped to 
// boundary condition objects), constructs a boundary cell map that provides
// access to metadata on the boundary cells (neighbor cell information, 
// face information, boundary conditions for faces).
boundary_cell_map_t* boundary_cell_map_from_mesh_and_bcs(mesh_t* mesh, string_ptr_unordered_map_t* bcs);

#endif
