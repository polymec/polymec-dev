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
