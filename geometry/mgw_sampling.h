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
                              double desired_spacing,
                              double ideal_energy_tolerance,
                              int* actual_num_points);

// Given a twice differentiable surface, this version computes the curvature 
// of the the surface and uses that to adjust the point spacing so that the 
// regions of higher curvature are better resolved.
point_t* adaptive_mgw_sampling(sp_func_t* C2_surface, 
                               bbox_t* bounding_box, 
                               int min_num_points,
                               double desired_spacing,
                               double ideal_energy_tolerance,
                               int* actual_num_points);

#endif

