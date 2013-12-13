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


#ifndef POLYMEC_FORCE_CVT_GEN_DIST_H
#define POLYMEC_FORCE_CVT_GEN_DIST_H

#include "geometry/cvt_gen_dist.h"

// This type represents a force exerted on a generator by neighboring 
// generators. Objects of this type are garbage-collected.
typedef struct cvt_gen_force_t cvt_gen_force_t;

// Creates a piecewise linear harmonic oscillator force model with the 
// given spring constant k and equilibrium displacement determined by the 
// density function. Here, F = k * (l - l0) for l < l0 and zero otherwise, 
// where l0 is the equilibrium length for a given region.
cvt_gen_force_t* linear_spring_force_new(double k);

// Creates a new force-balance-based CVT generator algorithm. The forces are 
// applied to the Delaunay triangulation of the generators to simplify 
// the topology. Arguments are:
// force - An object that computes the force between generators.
// tolerance - The threshold for the magnitude of the force at which the 
//             algorithm will terminate.
cvt_gen_dist_t* force_cvt_gen_dist_new(cvt_gen_force_t* force,
                                       double tolerance);

#endif

