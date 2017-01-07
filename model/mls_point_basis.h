// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_MLS_POINT_BASIS_H
#define POLYMEC_MLS_POINT_BASIS_H

#include "core/point_cloud.h"
#include "model/stencil.h"
#include "model/point_kernel.h"
#include "model/point_basis.h"

// Creates a moving-least-squares polynomial basis defined on the 
// given point cloud with neighborhoods given by the given stencil. Here,
// kernel_lengths is a field (array) assigning a characteristic extent, h, 
// to each point in the domain, and must have enough storage for ghost 
// points.
point_basis_t* mls_point_basis_new(int polynomial_order,
                                   point_kernel_t* kernel,
                                   point_cloud_t* domain,
                                   stencil_t* neighborhoods,
                                   real_t* kernel_lengths);

#endif
