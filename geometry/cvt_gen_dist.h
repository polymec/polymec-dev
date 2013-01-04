#ifndef POLYMEC_PROB_CVT_GEN_H
#define POLYMEC_PROB_CVT_GEN_H

#include "core/polymec.h"
#include "core/sp_func.h"

#ifdef __cplusplus
extern "C" {
#endif

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
void cvt_gen_dist_generate_random_points(long (*rng)(), sp_func_t* density, bbox_t* bounding_box, int num_points, point_t* points);

#ifdef __cplusplus
}
#endif

#endif

