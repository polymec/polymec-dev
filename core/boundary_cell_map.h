#include "core/mesh.h"
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

// This is an opaque type that represents a periodic boundary condition.
// Objects of this type are garbage-collected.
typedef struct periodic_bc_t periodic_bc_t;

// Creates a periodic boundary condition that identifies faces belonging 
// to tag1 with corresponding faces belonging to tag2. It is not necessary 
// that tag1[i] corresponds to tag2[i]--the geometric information about the 
// faces in the mesh will determine the correspondence.
periodic_bc_t* periodic_bc_new(const char* tag1, const char* tag2);

// Creates a periodic boundary condition with a function (and a context) 
// for generating periodic mappings between faces.
typedef int_int_unordered_map_t* (*periodic_bc_mapping_func)(void*, mesh_t*, char*, char*);
periodic_bc_t* periodic_bc_new_with_map_func(const char* tag1, const char* tag2, periodic_bc_mapping_func mapping_func, void* mapping_context);

// Returns true if the given periodic_bc object is actually a periodic_bc object,
// false if it has been improperly cast to one.
bool periodic_bc_is_valid(periodic_bc_t* bc);

// Returns true if the given pointer can be cast to a periodic_bc object, 
// false if not.
bool pointer_is_periodic_bc(void* ptr);

// Retrieves the two tags from the periodic boundary condition. tag1 and 
// tag2 will be pointed to internally-stored strings (which must not be deleted).
void periodic_bc_get_tags(periodic_bc_t* bc, char** tag1, char** tag2);

// Generates a face-to-face mapping for a periodic boundary condition.
int_int_unordered_map_t* periodic_bc_generate_map(periodic_bc_t* bc, mesh_t* mesh);

#ifdef __cplusplus
}
#endif

