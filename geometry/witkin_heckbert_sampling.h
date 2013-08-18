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
                                  double surface_diameter,
                                  int max_num_sample_points,
                                  point_t* initial_point,
                                  int* num_sample_points);

#endif

