// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_RECT_PRISM_SP_FUNC_H
#define POLYMEC_RECT_PRISM_SP_FUNC_H

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
sp_func_t* rect_prism_sp_func_new(point_t* x0, 
                                  real_t L1, real_t L2, real_t L3,
                                  real_t alpha, real_t beta, real_t gamma);

// Creates a rectangular prism identical to the given bounding box.
sp_func_t* rect_prism_sp_func_from_bbox(bbox_t* bounding_box);

#endif

