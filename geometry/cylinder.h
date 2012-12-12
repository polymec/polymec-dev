#ifndef POLYMEC_CYLINDER_H
#define POLYMEC_CYLINDER_H

#include "core/sp_func.h"

#ifdef __cplusplus
extern "C" {
#endif

// This signed distance function represents an infinite cylinder with a given 
// direction vector d, axial point x, and radius r. The orientation of the 
// normal vector (inward/outward) is also given.
sp_func_t* cylinder_new(vector_t* d, point_t* x, double r, normal_orient_t normal_orientation);

#ifdef __cplusplus
}
#endif

#endif

