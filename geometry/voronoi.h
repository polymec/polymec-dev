#ifndef ARBI_VORONOI_H
#define ARBI_VORONOI_H

#include "core/mesh.h"
#include "core/point.h"

#ifdef __cplusplus
extern "C" {
#endif

// This type represents a planar Voronoi graph.
typedef struct planar_voronoi_t planar_voronoi_t;

// This function creates a planar Voronoi graph given a set of Voronoi generators and an 
// implicit function F(x,y) which is zero at the domain boundary.
planar_voronoi_t* planar_voronoi_graph(point_t* points, int num_points, sp_func_t* F);

// This function creates a 3D Voronoi tessellation by extruding a planar Voronoi graph.
mesh_t* voronoi_extrusion(planar_voronoi_t* graph, int num_layers, double length);

// This function creates a Voronoi tessellation of the given points in 
// three-dimensional space, returning a fully-featured mesh.
mesh_t* voronoi_tessellation(point_t* points, int num_points, sp_func_t* F);

#ifdef __cplusplus
}
#endif

#endif

