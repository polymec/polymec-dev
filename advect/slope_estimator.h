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

#ifndef POLYMEC_SLOPE_ESTIMATOR_H
#define POLYMEC_SLOPE_ESTIMATOR_H

#include "core/point.h"

// This class computes a slope that is limited at discontinuities and is 
// equivalent to a centered difference where a quantity is continuous. It 
// constructs a parabolic least-squares fit to a solution in the vicinity of 
// a set of points, and uses that parabolic profile to determine the locations 
// of extrema. Objects of this type are garbage-collected.
typedef struct slope_estimator_t slope_estimator_t;

// Types of slope limiter functions.
typedef enum
{
  // These appear in order from least to most limited.
  SLOPE_LIMITER_BEAM_WARMING,  
  SLOPE_LIMITER_VAN_LEER,
  SLOPE_LIMITER_SUPERBEE,
  SLOPE_LIMITER_MINMOD,
  SLOPE_LIMITER_LAX_WENDROFF
} slope_limiter_t;

// Constructor for a slope estimator that uses the given type.
slope_estimator_t* slope_estimator_new(slope_limiter_t function);

// Computes the slope limiting factor on a given set of points, given a
// "center" point about which the fit is centered, a "forward" point that 
// determines the alignment of the least-squares parabolic fit for a 
// quantity, and the rest of the points.
// estimator        - The slope estimator object.
// center_point     - The point about which the parabolic least-squares fit is centered.
// center_value     - The value of the quantity at the center point.
// forward_point    - A point that can be used with the center point to construct 
//                    a finite difference estimate of the slope of the quantity.
// forward_value    - The value of the solution at the forward point.
// other_points     - An array of other points to use in the least-squares fit.
// other_values     - An array of solution values corresponding to the other points.
// num_other_points - The length of the other_points (and other_solutions) arrays.
double slope_estimator_value(slope_estimator_t* estimator, 
                             point_t* center_point,
                             double center_value,
                             point_t* forward_point,
                             double forward_value,
                             point_t* other_points,
                             double* other_values,
                             int num_other_points);

#endif
