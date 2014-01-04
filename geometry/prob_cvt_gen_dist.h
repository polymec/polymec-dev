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

#ifndef POLYMEC_PROB_CVT_GEN_DIST_H
#define POLYMEC_PROB_CVT_GEN_DIST_H

#include "geometry/cvt_gen_dist.h"

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
cvt_gen_dist_t* prob_cvt_gen_dist_new(int (*random_gen)(), 
                                      int num_samples, 
                                      real_t alpha, 
                                      real_t beta, 
                                      real_t min_dist,
                                      int max_iters);

#endif

