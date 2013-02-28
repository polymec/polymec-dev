#ifndef POLYMEC_CYLINDER_MAPPING_H
#define POLYMEC_CYLINDER_MAPPING_H

#include "core/sp_func.h"

#ifdef __cplusplus
extern "C" {
#endif

// This function maps a point xi within an indexed block to a point x in 
// 3D space. This mapping supports the construction of a block-structured 
// mesh consisting of five blocks that completely map a cylinder of given 
// radius and length. The block indices are as follows:
// 0: Inner block
// 1: -x block
// 2: +x block
// 3: -y block
// 4: +y block
// 5: -z block
// 6: +z block
// The cylinder is centered at the given axis that passes through the point 
// x0 and has radius r and length l. 
sp_func_t* cylinder_mapping_new(point_t* x0, vector_t* axis, double r, double l, int block_index);

#ifdef __cplusplus
}
#endif

#endif

