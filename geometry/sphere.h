#ifndef POLYMEC_SPHERE_H
#define POLYMEC_SPHERE_H

#include "core/sp_func.h"

#ifdef __cplusplus
extern "C" {
#endif

// This signed distance function represents a sphere centered at the given 
// point x with the given radius r. The normal orientation (outward/inward)
// is given by the third argument.
sp_func_t* sphere_new(point_t* x, double r, normal_orient_t normal_orientation);

#ifdef __cplusplus
}
#endif

#endif

