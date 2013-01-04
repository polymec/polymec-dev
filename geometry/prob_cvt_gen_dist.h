#ifndef POLYMEC_PROB_CVT_GEN_DIST_H
#define POLYMEC_PROB_CVT_GEN_DIST_H

#include "geometry/cvt_gen_dist.h"

#ifdef __cplusplus
extern "C" {
#endif

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
// min_dist    - Points that are closer to the boundary than this number are 
//               projected to the boundary.
// max_iters   - The maximum number of iterations to perform before quitting.
cvt_gen_dist_t* prob_cvt_gen_dist_new(long (*random_gen)(), 
                                      int num_samples, 
                                      double alpha, 
                                      double beta, 
                                      double min_dist,
                                      int max_iters);

#ifdef __cplusplus
}
#endif

#endif

