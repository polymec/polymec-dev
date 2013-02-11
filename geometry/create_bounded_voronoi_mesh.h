#ifndef POLYMEC_CREATE_BOUNDED_VORONOI_MESH_H
#define POLYMEC_CREATE_BOUNDED_VORONOI_MESH_H

#include "core/mesh.h"
#include "core/point.h"
#include "core/sp_func.h"

#ifdef __cplusplus
extern "C" {
#endif

// Creates a Voronoi tessellation of the given generator points in 
// three-dimensional space, bounded by a surface consisting of boundary 
// generators (cells whose generators lie on the boundary of the domain).
mesh_t* create_bounded_voronoi_mesh(point_t* generators, int num_generators,
                                    point_t* boundary_generators, int num_boundary_generators,
                                    point_t* ghost_generators, int num_ghost_generators);

// Creates a Voronoi tessellation of the given generator points in 
// three-dimensional space, bounded by an isosurface for which the given 
// "boundary" signed-distance function assumes the value zero.
//mesh_t* create_bounded_voronoi_mesh(point_t* generators, int num_generators,
//                                    point_t* ghost_generators, int num_ghost_generators,
//                                    sp_func_t* boundary);

#ifdef __cplusplus
}
#endif

#endif

