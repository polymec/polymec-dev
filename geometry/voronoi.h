#ifndef ARBI_VORONOI_H
#define ARBI_VORONOI_H

#include "core/mesh.h"
#include "core/point.h"
#include "core/sp_func.h"

#ifdef __cplusplus
extern "C" {
#endif

// This function creates a Voronoi tessellation of the given points in 
// three-dimensional space, returning a fully-featured mesh.
// The tessellation contains infinite cells at its boundary, and 
// these cells are tagged as "outer_cells". The edges bounding 
// these cells are semi-infinite, having one node each, and are 
// tagged as "outer_edges". These tags can be used to intersect 
// the mesh with a boundary function to form a closed domain.
mesh_t* voronoi_tessellation(point_t* points, int num_points, 
                             point_t* ghost_points, int num_ghost_points);

// This function creates a planar Voronoi tessellation, in which 
// a single layer of cells is bounded by two planes: one above, 
// and one below. The top plane is located at z = 0.5, and the 
// bottom one is located at z = -0.5.
mesh_t* voronoi_plane(point_t* points, int num_points, 
                      point_t* ghost_points, int num_ghost_points);

// This function intersects the given Voronoi tessellation (mesh) with the 
// given signed distance function representing a domain boundary, adding the 
// necessary faces and vertices. It is assumed that all generators fall within
// the domain.
void voronoi_intersect_with_boundary(mesh_t* mesh, sp_func_t* boundary);

// This function takes a planar Voronoi mesh and generates a 
// stack of these planes to form a new mesh. The bottom plane is 
// z = z1 and the top is z = z2.
mesh_t* voronoi_stack(mesh_t* plane, int num_planes, double z1, double z2);

// This function strips off "outer" cells and edges, leaving only the 
// interior set of cells.
void voronoi_prune(mesh_t* mesh);

#ifdef __cplusplus
}
#endif

#endif

