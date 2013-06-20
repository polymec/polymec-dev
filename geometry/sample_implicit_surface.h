#ifndef POLYMEC_SAMPLE_IMPLICIT_SURFACE_H
#define POLYMEC_SAMPLE_IMPLICIT_SURFACE_H

#include "core/sp_func.h"

// Given a desired density of sampling points per unit area (surface_density),
// this function computes and returns a set of points that sample the 
// given (closed) implicit surface. The surface diameter is an estimate of 
// the characteristic size of the surface, which is used to estimate a stopping 
// condition for the sampling. The algorithm terminates when it reaches the 
// desired surface density or the maximum number of sampling points. 
// The given initial point is used to seed the sampling process. 
// NOTE: if surface_density is NULL, a uniform density is assumed.
point_t* sample_implicit_surface(sp_func_t* surface, 
                                 sp_func_t* surface_density,
                                 double surface_diameter,
                                 int max_num_sample_points,
                                 point_t* initial_point,
                                 int* num_sample_points);

#endif

