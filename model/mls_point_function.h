// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_MLS_POINT_FUNCTION_H
#define POLYMEC_MLS_POINT_FUNCTION_H

#include "core/point_cloud.h"
#include "model/stencil.h"
#include "model/point_kernel.h"
#include "model/point_function.h"

// Creates a moving-least-squares polynomial shape function defined on the 
// given point cloud with neighborhoods given by the given stencil. Here,
// smoothing_lengths is a field (array) assigning a characteristic extent, h, 
// to each point in the domain, and must have enough storage for ghost 
// points.
point_function_t* mls_point_function_new(int polynomial_order,
                                         point_kernel_t* kernel,
                                         point_cloud_t* domain,
                                         stencil_t* neighborhoods,
                                         real_t* kernel_lengths);

#endif
