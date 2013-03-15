#ifndef POLYMEC_HEXAHEDRON_H
#define POLYMEC_HEXAHEDRON_H

#include "core/point.h"

#ifdef __cplusplus
extern "C" {
#endif

// A hexahedron is a deformable region that is topologically hexahedral. Its 
// geometry is defined by a set of nodes in 3-dimensional space.
// Objects of this type are garbage-collected.
//
// For nodes, we use the reference cell below, which is typical of 3D
// finite element schemes:
//             
//     7o----6o      z^  y
//     /|    /|       | /
//   4o----5o |       |/   x
//    |3o---|2o       +---->
//    |/    |/       
//   0o----1o      
//
// Hexahedra have a geometric "order" that determines the number of nodes 
// along their edges and within their interior. The nodes are placed on a 
// cubic lattice within the logical space of the hexahedron, and there are 
// (s+1) nodes in each of the 3 directions, where s is the geometric order.
// 
// The faces are numbered 0-5, with 0-1 being the (-/+) x faces,
// 2-3 the (-/+) y faces, and 4-5 the (-/+) z faces.
//
// Block faces are indexed as follows:
//  0 - negative 'x' face
//  1 - positive 'x' face
//  2 - negative 'y' face
//  3 - positive 'y' face
//  4 - negative 'z' face
//  5 - positive 'z' face
// Note that the coordinates are in quotes because 'x', 'y', and 'z' are 
// topological, not geometric coordinates. In other words, 'x' may not be 
// aligned with x in physical space, etc.
typedef struct hexahedron_t hexahedron_t;

// Creates a hexahedron of the given order.
hexahedron_t* hexahedron_new(int order);

// Returns the geometric order of the hexahedron.
int hexahedron_order(hexahedron_t* hex);

// Returns the number of points that define the geometry/topology of the 
// hexahedron. This is related to the order of the hexahedron, but provides a 
// shortcut.
int hexahedron_num_points(hexahedron_t* hex);

// Retrieves the logical coordinates of the points within the hexahedron.
void hexahedron_get_points(hexahedron_t* hex, point_t* points);

// Maps a point xi within the logical space of the hexahedron to its 
// equivalent point x in physical space, using the isoparametric mapping 
// for the hexahedron defined by the given points.
void hexahedron_map(hexahedron_t* hex, point_t* points, point_t* xi, point_t* x);

#ifdef __cplusplus
}
#endif

#endif

