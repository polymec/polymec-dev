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

#ifndef POLYMEC_MGW_SAMPLING_H
#define POLYMEC_MGW_SAMPLING_H

#include "core/sp_func.h"

// These functions compute and returns a set of points that sample the given 
// (closed) implicit surface using the method of Meyer, Georgel, and Whitaker 
// (Shape Modeling International, 2005). For each function, the given bounding 
// box is used to initialize initial (random) point positions and should have 
// dimensions roughly matching that of the surface. The implicit function 
// representing the surface must have at least one spatial derivative.

// This version produces the given number of uniform sample points on the 
// given surface. The algorithm guarantees a minimum number of sample points 
// and a desired spacing for the points, adding or deleting them where 
// necessary. A relative tolerance is used to gauge whether the energy of the 
// points is close enough to the "ideal energy" that produces a hexagonal 
// close packing of points on the surface.
point_t* uniform_mgw_sampling(sp_func_t* C1_surface, 
                              bbox_t* bounding_box, 
                              int min_num_points,
                              real_t desired_spacing,
                              real_t ideal_energy_tolerance,
                              int* actual_num_points);

// Given a twice differentiable surface, this version computes the curvature 
// of the the surface and uses that to adjust the point spacing so that the 
// regions of higher curvature are better resolved.
point_t* adaptive_mgw_sampling(sp_func_t* C2_surface, 
                               bbox_t* bounding_box, 
                               int min_num_points,
                               real_t desired_spacing,
                               real_t ideal_energy_tolerance,
                               int* actual_num_points);

#endif

