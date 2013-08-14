#ifndef POLYMEC_CREATE_VORONOI_MESH_H
#define POLYMEC_CREATE_VORONOI_MESH_H

#include "core/mesh.h"
#include "core/point.h"
#include "core/slist.h"

// This function creates a Voronoi tessellation of the given points in 
// three-dimensional space. The tessellation contains only bounded cells.
// This means that generators corresponding to unbounded cells do not appear 
// in the tessellation, so it's important to include any generators that 
// produce a non-zero set of bounded Voronoi cells. Optionally, a linked list 
// may be provided that will store the indices of any cells that were deleted to construct a 
// completely bounded mesh. If such information is not desired, the final 
// argument can be set to NULL.
mesh_t* create_voronoi_mesh(point_t* generators, int num_generators, 
                            point_t* ghost_generators, int num_ghost_generators,
                            int_slist_t* deleted_cells);

#endif

