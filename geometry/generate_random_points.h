// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_GENERATE_RANDOM_POINTS_H
#define POLYMEC_GENERATE_RANDOM_POINTS_H

#include "core/sp_func.h"

/// \addtogroup geometry geometry
///@{

/// This function generates the given number of points within the
/// given bounding box, from the given probability density function. The
/// given random number generator is used.
void generate_random_points(rng_t* rng, sp_func_t* density, bbox_t* bounding_box, size_t num_points, point_t* points);

///@}

#endif

