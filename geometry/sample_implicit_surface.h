#ifndef POLYMEC_SAMPLE_IMPLICIT_SURFACE_H
#define POLYMEC_SAMPLE_IMPLICIT_SURFACE_H

#include "core/sp_func.h"

// Given a desired density of sampling points per unit area (surface_density),
// this function computes and returns a set of points that sample the 
// given (closed) implicit surface. The given bounding box must bound the 
// surface and is used to generate the random starting point.
point_t* sample_implicit_surface(sp_func_t* surface, 
                                 sp_func_t* surface_density,
                                 bbox_t* bounding_box,
                                 int* num_sample_points);

#endif

