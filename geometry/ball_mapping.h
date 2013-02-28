#ifndef POLYMEC_BALL_MAPPING_H
#define POLYMEC_BALL_MAPPING_H

#include "core/sp_func.h"

#ifdef __cplusplus
extern "C" {
#endif

// This function maps a point xi within an indexed block to a point x in 
// 3D space. This mapping supports the construction of a block-structured 
// mesh consisting of seven blocks that completely map a ball of given 
// radius. The block indices are as follows:
// 0: Inner block
// 1: -x block
// 2: +x block
// 3: -y block
// 4: +y block
// 5: -z block
// 6: +z block
// The ball is centered at x0 and has radius r. The length on a side of the 
// inner block is l, which must be less than 2*r.
sp_func_t* ball_mapping_new(point_t* x0, double r, double l, int block_index);

#ifdef __cplusplus
}
#endif

#endif

