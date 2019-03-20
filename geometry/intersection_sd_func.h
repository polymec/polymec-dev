// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_INTERSECTION_SD_FUNC_H
#define POLYMEC_INTERSECTION_SD_FUNC_H

#include "geometry/sd_func.h"

/// \addtogroup geometry geometry
///@{

/// This signed distance function represents the intersection of the set of
/// given surfaces represented by signed distance functions.
/// \relates sd_func
sd_func_t* intersection_sd_func_new(sd_func_t** surfaces, int num_surfaces);

/// This is a shorthand function that creates the intersection of two
/// surfaces.
/// \relates sd_func
sd_func_t* intersection_sd_func_new2(sd_func_t* surface1, sd_func_t* surface2);

///@}

#endif

