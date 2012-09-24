#ifndef ARBI_CREATE_CUBIC_LATTICE_MESH_H
#define ARBI_CREATE_CUBIC_LATTICE_MESH_H

#include "core/mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

// This function creates and returns a mesh for a cubic lattice of 
// nx x ny x nz cells. The mesh spans the interval [0,1] x [0,1] x [0,1].
mesh_t* create_cubic_lattice_mesh(int nx, int ny, int nz);

#ifdef __cplusplus
}
#endif

#endif

