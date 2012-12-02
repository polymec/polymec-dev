#ifndef POLYMEC_BOUNDARY_CELL_MAP_H
#define POLYMEC_BOUNDARY_CELL_MAP_H

#include "core/mesh.h"
#include "core/periodic_bc.h"
#include "core/unordered_map.h"

#ifdef __cplusplus
extern "C" {
#endif

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
  // entries identifying the cells on the opposite side of the face.
  // opp_cells[i] is NULL if the ith boundary face does not have a periodic 
  // boundary condition.
  cell_t** opp_cells;
} boundary_cell_t;

// Define the mapping type.
DEFINE_UNORDERED_MAP(boundary_cell_map, int, boundary_cell_t*, int_hash, int_equals)

// Given a mesh and a table of boundary conditions (mesh tags mapped to 
// boundary condition objects), constructs a boundary cell map that provides
// access to metadata on the boundary cells (neighbor cell information, 
// face information, boundary conditions for faces).
boundary_cell_map_t* boundary_cell_map_from_mesh_and_bcs(mesh_t* mesh, str_ptr_unordered_map_t* bcs);

#ifdef __cplusplus
}
#endif

#endif
