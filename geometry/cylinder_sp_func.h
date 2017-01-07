// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_CYLINDER_SP_FUNC_H
#define POLYMEC_CYLINDER_SP_FUNC_H

#include "core/sp_func.h"

// This signed distance function represents an infinite cylinder with a given 
// axial point x and radius r. The orientation of the normal vector 
// (inward/outward) is also given.
sp_func_t* cylinder_sp_func_new(point_t* x, real_t r, normal_orient_t normal_orientation);

#endif

