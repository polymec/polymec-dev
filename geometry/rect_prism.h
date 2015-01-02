// Copyright (c) 2012-2015, Jeffrey N. Johnson
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
                          real_t L1, real_t L2, real_t L3,
                          real_t alpha, real_t beta, real_t gamma);

// Creates a rectangular prism identical to the given bounding box.
sp_func_t* rect_prism_from_bbox(bbox_t* bounding_box);

#endif

