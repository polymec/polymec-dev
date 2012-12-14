#ifndef POLYMEC_PROB_CVT_GEN_H
#define POLYMEC_PROB_CVT_GEN_H

#include "core/polymec.h"
#include "core/sp_func.h"

#ifdef __cplusplus
extern "C" {
#endif

// This class creates a set of Voronoi generators that are properly placed 
// to create a centroidal Voronoi tessellation, using a probabilistic algorithm 
// described in Du and Gunzburger, Parallel Comp. 28 (2002) pp. 1477-1500.
// This is a fast approximate CVT generation method that does NOT require the 
// Voronoi tessellation of the generators to iterate the points, and is thus 
// easily parallelized. Objects of this type are garbage-collected.
typedef struct prob_cvt_gen_t prob_cvt_gen_t;

// This class provides an interface for terminating the iteration of a 
// probabilistic CVT generator algorithm. Objects of this type are garbage-
// collected.
typedef struct prob_cvt_gen_term_t prob_cvt_gen_term_t;

// Creates a new probabilistic CVT generator algorithm with the following 
// parameters:
// random_gen  - A function that returns a random integer between 0 and RAND_MAX.
// num_samples - The number of sampling points used in an iteration to correct
//               each generator position.
// alpha, beta - Coefficients for the generator position correction algorithm,
//               which is: 
// 
//             (alpha*ji + beta)*zi + ((1-alpha)*ji + (1-beta))*ui
//     zi <--- ---------------------------------------------------
//                                ji + 1
//
//               where zi is the position of the ith generator, ui is the 
//               average position of the set of sample points within the 
//               voronoi region of generator i, and ji is an iteration index
//               that forces the algorithm to converge by weighting either 
//               zi or ui more and more as the iteration number increases.
prob_cvt_gen_t* prob_cvt_gen_new(long (*random_gen)(), int num_samples, double alpha, double beta);

// Given an initial set of generator points, move them around according to 
// the algorithm described above until termination critierion (provided by 
// termination) is achieved. The density function is used to create the 
// sampling points. Boundary is a signed distance function which is negative 
// inside the boundary, zero on the surface, and positive outside.
// NOTE: boundary can be NULL, but bounding_box must be given. 
// NOTE: If boundary is NULL, the bounding box is assumed to describe a 
// NOTE: rectangular domain. 
void prob_cvt_gen_iterate(prob_cvt_gen_t* prob, 
                          sp_func_t* density,
                          sp_func_t* boundary,
                          bbox_t* bounding_box,
                          prob_cvt_gen_term_t* termination,
                          point_t* points, 
                          int num_points);

// Creates an object that terminates a probabilistic CVT iteration after 
// the given number of iterations.
prob_cvt_gen_term_t* terminate_prob_cvt_at_iter(int max_iter);

#ifdef __cplusplus
}
#endif

#endif

