// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_SCALED_SD_FUNC_H
#define POLYMEC_SCALED_SD_FUNC_H

#include "geometry/sd_func.h"

/// \addtogroup geometry geometry
///@{

/// This signed distance function takes another such function and scales it
/// by a given factor.
/// \relates sd_func
sd_func_t* scaled_sd_func_new(sd_func_t* func, real_t scale_factor);

///@}

#endif

