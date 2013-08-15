#ifndef POLYMEC_CREATE_BOUNDARY_GENERATORS_H
#define POLYMEC_CREATE_BOUNDARY_GENERATORS_H

#include "core/point.h"
#include "core/array.h"

// This function creates a set of generators that can be used as stationary 
// generators in a Centroidal Voronoi Tessellation (CVT) algorithm.
void create_boundary_generators(ptr_array_t* surface_points, 
                                ptr_array_t* surface_normals, 
                                ptr_array_t* surface_tags,
                                point_t** boundary_generators,
                                int* num_boundary_generators,
                                char*** tag_names,
                                int_array_t*** tags,
                                int* num_tags);

#endif

