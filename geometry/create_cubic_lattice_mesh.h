#ifndef POLYMEC_CREATE_CUBIC_LATTICE_MESH_H
#define POLYMEC_CREATE_CUBIC_LATTICE_MESH_H

#include "core/mesh.h"

// This function creates and returns a mesh for a cubic lattice of 
// nx x ny x nz cells. The mesh spans the rectangular region of space 
// defined by the bounding box.
mesh_t* create_cubic_lattice_mesh_with_bbox(int nx, int ny, int nz, bbox_t* bbox);

// This function creates and returns a mesh for a cubic lattice of 
// nx x ny x nz cells. The mesh spans the interval [0,1] x [0,1] x [0,1].
mesh_t* create_cubic_lattice_mesh(int nx, int ny, int nz);

// This function tags the faces of a newly-created cubic lattice mesh 
// for convenient boundary condition assignments.
void tag_cubic_lattice_mesh_faces(mesh_t* mesh, int nx, int ny, int nz,
                                  const char* x1_tag, 
                                  const char* x2_tag, 
                                  const char* y1_tag,
                                  const char* y2_tag,
                                  const char* z1_tag,
                                  const char* z2_tag);

#endif

