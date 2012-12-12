#ifndef POLYMEC_VORONOI_H
#define POLYMEC_VORONOI_H

#include "core/mesh.h"
#include "core/faceted_surface.h"
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
mesh_t* unbounded_voronoi(point_t* generators, int num_generators, 
                          point_t* ghost_generators, int num_ghost_generators);

// Creates a Voronoi tessellation of the given generator points in 
// three-dimensional space, bounded by an isosurface for which the given 
// "boundary" signed-distance function assumes the value zero.
mesh_t* bounded_voronoi(point_t* generators, int num_generators,
                        point_t* ghost_generators, int num_ghost_generators,
                        sp_func_t* boundary);

#if 0
// This function intersects the given Voronoi tessellation (mesh) with the 
// given signed distance function representing a domain boundary, returning 
// a faceted surface representing the intersection.
faceted_surface_t* voronoi_intersect_with_boundary(mesh_t* mesh, sp_func_t* boundary);

// This function strips off the infinite "outer" cells and edges generated 
// by voronoi_tessellation (those tagged as "outer_cells" and "outer_edges").
void voronoi_prune(mesh_t* mesh);

// This function creates a Voronoi tessellation that is contained within 
// the given faceted surface.
mesh_t* voronoi_tessellation_within_surface(point_t* points, int num_points,
                                            point_t* ghost_points, int num_ghost_points,
                                            faceted_surface_t* surface);

// This function creates a planar Voronoi tessellation, in which 
// a single layer of cells is bounded by two planes: one above, 
// and one below. The top plane is located at z = 0.5, and the 
// bottom one is located at z = -0.5.
mesh_t* voronoi_plane(point_t* points, int num_points, 
                      point_t* ghost_points, int num_ghost_points);

// This function takes a planar Voronoi mesh and generates a 
// stack of these planes to form a new mesh. The bottom plane is 
// z = z1 and the top is z = z2.
mesh_t* voronoi_stack(mesh_t* plane, int num_planes, double z1, double z2);
#endif

#ifdef __cplusplus
}
#endif

#endif

