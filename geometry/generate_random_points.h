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

#ifndef POLYMEC_GENERATE_RANDOM_POINTS_H
#define POLYMEC_GENERATE_RANDOM_POINTS_H

#include "core/sp_func.h"

// This function generates the given number of points within the 
// given bounding box, from the given probability density function. The 
// given random number generator is used.
void generate_random_points(long (*rng)(), sp_func_t* density, bbox_t* bounding_box, int num_points, point_t* points);

#endif

