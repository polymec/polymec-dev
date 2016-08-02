// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_SHEPARD_POINT_BASIS_H
#define POLYMEC_SHEPARD_POINT_BASIS_H

#include "core/point_cloud.h"
#include "model/stencil.h"
#include "model/point_kernel.h"
#include "model/point_basis.h"

// Creates a Shepard basis function defined on the given point cloud with 
// neighborhoods given by the given stencil. Here, kernel_lengths is a 
// field (array) assigning a characteristic extent, h, to each point in the 
// domain, and must have enough storage for ghost points.
point_basis_t* shepard_point_basis_new(point_kernel_t* kernel,
                                       point_cloud_t* domain,
                                       stencil_t* neighborhoods,
                                       real_t* kernel_lengths);

#endif
