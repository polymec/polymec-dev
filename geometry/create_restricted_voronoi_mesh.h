#ifndef POLYMEC_CREATE_RESTRICTED_VORONOI_MESH_H
#define POLYMEC_CREATE_RESTRICTED_VORONOI_MESH_H

#include "core/mesh.h"
#include "core/surface_mesh.h"
#include "core/point_set.h"

#ifdef __cplusplus
extern "C" {
#endif

// Creates a Voronoi tessellation of the given generator points in 
// three-dimensional space, restricted by the given surface mesh.
// The generators are given in a point_set (K-d tree).
// generators (cells whose generators lie on the boundary of the domain).
mesh_t* create_restricted_voronoi_mesh(point_set_t* generators,
                                       point_set_t* ghost_generators,
                                       surface_mesh_t* surface_mesh);

#ifdef __cplusplus
}
#endif

#endif

