#ifndef POLYMEC_GENERATE_RANDOM_POINTS_H
#define POLYMEC_GENERATE_RANDOM_POINTS_H

#include "core/sp_func.h"

#ifdef __cplusplus
extern "C" {
#endif

// This function generates the given number of points within the 
// given bounding box, from the given probability density function. The 
// given random number generator is used.
void generate_random_points(long (*rng)(), sp_func_t* density, bbox_t* bounding_box, int num_points, point_t* points);

#ifdef __cplusplus
}
#endif

#endif

