#ifndef POLYMEC_CREATE_CUBIC_LATTICE_MESH_H
#define POLYMEC_CREATE_CUBIC_LATTICE_MESH_H

#include "core/mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

// This function creates and returns a mesh for a cubic lattice of 
// nx x ny x nz cells. The mesh spans the rectangular region of space 
// defined by the bounding box.
mesh_t* create_cubic_lattice_mesh_with_bbox(int nx, int ny, int nz, bbox_t* bbox);

// This function creates and returns a mesh for a cubic lattice of 
// nx x ny x nz cells. The mesh spans the interval [0,1] x [0,1] x [0,1].
mesh_t* create_cubic_lattice_mesh(int nx, int ny, int nz);

#ifdef __cplusplus
}
#endif

#endif

