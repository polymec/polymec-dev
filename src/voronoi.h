#ifndef ARBI_VORONOI_H
#define ARBI_VORONOI_H

#include "mesh.h"
#include "point.h"

#ifdef __cplusplus
extern "C" {
#endif

// This function creates a Voronoi tessellation of the given points in 
// three-dimensional space, returning a fully-featured mesh.
mesh_t* voronoi_tessellation(point_t* points, int num_points, bbox_t* bounding_box);

// This function creates a planar Voronoi graph using Fortune's algorithm. 
// Z coordinates are set to zero,and the cells have unit length in the z direction.
mesh_t* planar_voronoi_tessellation(point_t* points, int num_points, bbox_t* bounding_box);

// This function produces a 3D mesh by linearly extruding the given planar graph 
// in the z direction.
mesh_t* extrusion(mesh_t* planar_mesh, int num_extruded_cells, double length);

#ifdef __cplusplus
}
#endif

#endif

