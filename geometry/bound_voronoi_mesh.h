#ifndef POLYMEC_BOUND_VORONOI_MESH_H
#define POLYMEC_BOUND_VORONOI_MESH_H

#include "core/mesh.h"
#include "core/sp_func.h"
#include "core/mesh_diff.h"

// This takes an unbounded Voronoi tessellation/mesh and generates a set of 
// changes that bound it with the boundary represented by the given implicit 
// function. These changes are stored in the mesh_diff object that is returned.
mesh_diff_t* bound_voronoi_mesh(mesh_t* mesh, sp_func_t* boundary);

#endif

