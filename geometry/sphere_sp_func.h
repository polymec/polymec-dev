// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_SPHERE_SP_FUNC_H
#define POLYMEC_SPHERE_SP_FUNC_H

#include "core/sp_func.h"

// This signed distance function represents a sphere centered at the given 
// point x with the given radius r. The normal orientation (outward/inward)
// is given by the third argument.
sp_func_t* sphere_sp_func_new(point_t* x, real_t r, normal_orient_t normal_orientation);

#endif

