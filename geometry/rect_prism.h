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

#ifndef POLYMEC_RECT_PRISM_H
#define POLYMEC_RECT_PRISM_H

#include "core/sp_func.h"

// This signed distance function represents a rectangular prism that is 
// arbitrarily rotated in space. Arguments:
// x0 - The center point of the rectangular prism.
// L1, L2, L3 - Spatial extents of the primary, secondary, and tertiary axes 
//              of the rectangular prism.
// alpha, beta, gamma - Euler angles identifying the axes X, Y, Z of the 
//                      prism with respect to the usual Cartesian (x, y, z) 
//                      coordinate axes:
//              alpha - A rotation about the z axis.
//               beta - A rotation about the N axis (line of nodes).
//              gamma - A rotation about the Z axis.
sp_func_t* rect_prism_new(point_t* x0, 
                          double L1, double L2, double L3,
                          double alpha, double beta, double gamma);

// Creates a rectangular prism identical to the given bounding box.
sp_func_t* rect_prism_new_from_bbox(bbox_t* bounding_box);

#endif

