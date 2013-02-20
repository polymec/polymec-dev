#ifndef POLYMEC_CROP_VORONOI_MESH_H
#define POLYMEC_CROP_VORONOI_MESH_H

#include "core/mesh.h"
#include "core/surface_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

// Given a Voronoi tessellation/mesh, crop it with the given surface 
// mesh by chopping the open cells at the boundary.
void crop_voronoi_mesh(mesh_t* mesh, surface_mesh_t* surface_mesh);

#ifdef __cplusplus
}
#endif

#endif

