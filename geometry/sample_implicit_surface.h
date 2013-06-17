#ifndef POLYMEC_SAMPLE_IMPLICIT_SURFACE_H
#define POLYMEC_SAMPLE_IMPLICIT_SURFACE_H

#include "core/point.h"

// Given a desired density of sampling points per unit area (surface_density),
// this function computes and returns a set of points that sample the 
// given (closed) implicit surface.
point_t* sample_implicit_surface(sp_func_t* surface_density, 
                                 sp_func_t* surface,
                                 int* num_sample_points);

#endif

