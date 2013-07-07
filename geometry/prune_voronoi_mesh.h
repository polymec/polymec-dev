#ifndef POLYMEC_PRUNE_VORONOI_MESH_H
#define POLYMEC_PRUNE_VORONOI_MESH_H

#include "core/mesh.h"

// This takes an unbounded Voronoi tessellation/mesh and prunes the "outer" 
// (semi-infinite) cells and edges from it. This is not a reversible 
// operation.
void prune_voronoi_mesh(mesh_t* mesh);

#endif

