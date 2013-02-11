#ifndef POLYMEC_CREATE_UNBOUNDED_VORONOI_MESH_H
#define POLYMEC_CREATE_UNBOUNDED_VORONOI_MESH_H

#include "core/mesh.h"
#include "core/point.h"
#include "core/sp_func.h"

#ifdef __cplusplus
extern "C" {
#endif

// This function creates a Voronoi tessellation of the given points in 
// three-dimensional space. The tessellation contains infinite cells at its 
// boundary, and these cells are tagged as "outer_cells". The edges bounding 
// these cells are semi-infinite, having one node each, and are tagged as 
// "outer_edges". These tags can be used to intersect the mesh with a 
// boundary function to form a closed domain.
mesh_t* create_unbounded_voronoi_mesh(point_t* generators, int num_generators, 
                                      point_t* ghost_generators, int num_ghost_generators);

#ifdef __cplusplus
}
#endif

#endif

