#ifndef POLYMEC_CREATE_CVT_WITH_LLOYD_ITERATION_H
#define POLYMEC_CREATE_CVT_WITH_LLOYD_ITERATION_H

#include "geometry/create_cvt.h"

// Creates a Centroidal Voronoi Tessellation using Lloyd iteration with 
// the specified number of iterations.
mesh_t* create_cvt_with_lloyd_iteration(point_t* stationary_generators, int num_stationary_generators, 
                                        point_t* mobile_generators, int num_mobile_generators,
                                        int num_iterations);

#endif

