// Copyright 2012-2013 Jeffrey Johnson.
// 
// This file is part of Polymec, and is licensed under the Apache License, 
// Version 2.0 (the "License"); you may not use this file except in 
// compliance with the License. You may may find the text of the license in 
// the LICENSE file at the top-level source directory, or obtain a copy of 
// it at
// 
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef POLYMEC_BOUNDARY_CELL_MAP_H
#define POLYMEC_BOUNDARY_CELL_MAP_H

#include "core/mesh.h"
#include "core/periodic_bc.h"
#include "core/unordered_map.h"

// This data structure provides metadata for cells that touch the boundary 
// of a domain.
typedef struct
{
  // Cells neighboring this boundary cell.
  int num_neighbor_cells;
  int* neighbor_cells;

  // Boundary faces attached to this cell.
  int num_boundary_faces;
  int* boundary_faces;

  // bc_for_face[i] stores a pointer to the boundary condition for the 
  // ith boundary face (boundary_faces[i]) of this cell.
  // NOTE: If the ith boundary face has a periodic boundary condition,
  // NOTE: boundary_faces[i] will be NULL.
  void** bc_for_face;

  // For faces with periodic boundary conditions, this array has non-NULL
  // entries identifying the face on the opposite side of the ith face.
  // opp_faces[i] is NULL if the ith boundary face does not have a periodic 
  // boundary condition. 
  face_t** opp_faces;
} boundary_cell_t;

// Define the mapping type.
DEFINE_UNORDERED_MAP(boundary_cell_map, int, boundary_cell_t*, int_hash, int_equals)

// Given a mesh and a table of boundary conditions (mesh tags mapped to 
// boundary condition objects), constructs a boundary cell map that provides
// access to metadata on the boundary cells (neighbor cell information, 
// face information, boundary conditions for faces).
boundary_cell_map_t* boundary_cell_map_from_mesh_and_bcs(mesh_t* mesh, string_ptr_unordered_map_t* bcs);

#endif
