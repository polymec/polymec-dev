// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_SPHERE_SD_FUNC_H
#define POLYMEC_SPHERE_SD_FUNC_H

#include "geometry/sd_func.h"

/// \addtogroup geometry geometry
///@{

/// This signed distance function represents a sphere centered at the given
/// point x with the given radius r. The normal orientation (outward/inward)
/// is given by the third argument.
/// \relates sd_func
sd_func_t* sphere_sd_func_new(point_t* x, real_t r, normal_orient_t normal_orientation);

///@}

#endif

