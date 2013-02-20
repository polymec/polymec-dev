#ifndef POLYMEC_CREATE_BLOCK_MESH_H
#define POLYMEC_CREATE_BLOCK_MESH_H

#include "core/mesh.h"
#include "geometry/block.h"

#ifdef __cplusplus
extern "C" {
#endif

// Creates an unstructured mesh whose cells are hexahedrons contained 
// within logically hexahedral blocks given by a block assembly. The blocks 
// are enumerated within the assembly, and their corresponding resolutions 
// are specified by the array of block grids.
mesh_t* create_block_mesh(block_assembly_t* assembly, block_grid_t* resolutions);

#ifdef __cplusplus
}
#endif

#endif

