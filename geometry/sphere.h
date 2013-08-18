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

#ifndef POLYMEC_SPHERE_H
#define POLYMEC_SPHERE_H

#include "core/sp_func.h"

// This signed distance function represents a sphere centered at the given 
// point x with the given radius r. The normal orientation (outward/inward)
// is given by the third argument.
sp_func_t* sphere_new(point_t* x, double r, normal_orient_t normal_orientation);

#endif

