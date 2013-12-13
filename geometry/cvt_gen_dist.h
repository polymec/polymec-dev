// Copyright (c) 2012-2013, Jeffrey N. Johnson
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this 
// list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice, 
// this list of conditions and the following disclaimer in the documentation 
// and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#ifndef POLYMEC_PROB_CVT_GEN_H
#define POLYMEC_PROB_CVT_GEN_H

#include "core/polymec.h"
#include "core/sp_func.h"

// This base class provides an interface for algorithms that create a set of 
// Voronoi generators for a centroidal Voronoi tessellation. Objects of this 
// type are garbage-collected.
typedef struct cvt_gen_dist_t cvt_gen_dist_t;

// A function pointer type for performing the iteration.
typedef void (*cvt_gen_dist_iterate_func)(void*, sp_func_t*, sp_func_t*, 
                                          bbox_t*, point_t*, int);

// A destructor for any given context object.
typedef void (*cvt_gen_dist_dtor)(void*);

// This virtual table must be implemented by any CVT generator distribution
// algorithm.
typedef struct 
{
  cvt_gen_dist_iterate_func     iterate;
  cvt_gen_dist_dtor             dtor;
} cvt_gen_dist_vtable;

// Creates a new CVT generator distribution algorithm.
cvt_gen_dist_t* cvt_gen_dist_new(const char* name, void* context, cvt_gen_dist_vtable vtable);

// Returns the name of the CVT generator distribution algorithm.
const char* cvt_gen_dist_name(cvt_gen_dist_t* dist);

// Given an initial set of generator points, move them around according to 
// the designated algorithm until some termination critierion is achieved. 
// The density function is a relative measure of the number of generators per 
// unit volume. The boundary signed distance function is negative inside the 
// boundary, zero on the surface, and positive outside. The points that 
// fall on the boundary are moved to the end of the list of points, and the 
// number of boundary points is stored in num_boundary_points.
// NOTE: boundary can be NULL, but bounding_box must be given. 
// NOTE: If boundary is NULL, the bounding box is assumed to describe a 
// NOTE: rectangular domain. 
void cvt_gen_dist_iterate(cvt_gen_dist_t* dist, 
                          sp_func_t* density,
                          sp_func_t* boundary,
                          bbox_t* bounding_box,
                          point_t* points, 
                          int num_points,
                          int* num_boundary_points);

// This helper function generates the given number of points within the 
// given bounding box, from the given probability density function. The 
// given random number generator is used.
void cvt_gen_dist_generate_random_points(int (*rng)(), sp_func_t* density, bbox_t* bounding_box, int num_points, point_t* points);

#endif

