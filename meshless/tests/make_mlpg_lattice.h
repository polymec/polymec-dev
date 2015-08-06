// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/point_cloud.h"
#include "model/stencil.h"

// Creates a point cloud, a set of extents, and (optionally) a stencil that 
// can be used with MLPG-based methods.
void make_mlpg_lattice(int nx, int ny, int nz, real_t R_over_dx,
                       point_cloud_t** domain,
                       real_t** extents,
                       stencil_t** neighborhoods);

