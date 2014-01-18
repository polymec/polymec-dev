// Copyright (c) 2012-2014, Jeffrey N. Johnson
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

#ifndef POLYMEC_WITKIN_HECKBERT_SAMPLING_H
#define POLYMEC_WITKIN_HECKBERT_SAMPLING_H

#include "core/sp_func.h"

// Given a desired density of sampling points per unit area (surface_density),
// this function computes and returns a set of points that sample the 
// given (closed) implicit surface using the method of Witkin and Heckbert 
// (Proceedings of the 21st annual conference on Computer graphics and 
// interactive techniques, pp. 269-277, ACM Press, 1994). The surface diameter 
// is an estimate of the characteristic size of the surface, which is used to 
// estimate a stopping condition for the sampling. The algorithm terminates 
// when it reaches the desired surface density or the maximum number of 
// sampling points. The given initial point is used to seed the sampling 
// process. 
// NOTE: if surface_density is NULL, a uniform density is assumed.
point_t* witkin_heckbert_sampling(sp_func_t* surface, 
                                  sp_func_t* surface_density,
                                  real_t surface_diameter,
                                  int max_num_sample_points,
                                  point_t* initial_point,
                                  int* num_sample_points);

#endif

