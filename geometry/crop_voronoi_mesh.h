#ifndef POLYMEC_CROP_VORONOI_MESH_H
#define POLYMEC_CROP_VORONOI_MESH_H

#include "core/mesh.h"
#include "core/surface_mesh.h"
#include "core/mesh_diff.h"

// Given a Voronoi tessellation/mesh, generate a mesh_diff containing a set 
// of changes to a mesh that crop it with the given surface mesh.
mesh_diff_t* crop_voronoi_mesh(mesh_t* mesh, surface_mesh_t* surface_mesh);

#endif

