// Copyright 2012-2013 Jeffrey Johnson.
// 
// This file is part of Polymec, and is licensed under the Apache License, 
// Version 2.0 (the "License"); you may not use this file except in 
// compliance with the License. You may may find the text of the license in 
// the LICENSE file at the top-level source directory, or obtain a copy of 
// it at
// 
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

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

