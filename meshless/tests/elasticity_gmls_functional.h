// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ELASTICITY_GMLS_FUNCTIONAL_H
#define ELASTICITY_GMLS_FUNCTIONAL_H

#include "core/point_cloud.h"
#include "meshless/gmls_functional.h"

// Constructs a GMLS functional for solving the equations of linear elasticity 
// on the interior of a domain, using a polynomial basis of the given degree.
// Material parameters are Young's modulus E and Poisson's ratio nu.
gmls_functional_t* elasticity_gmls_functional_new(real_t E, real_t nu,
                                                  int degree,
                                                  point_cloud_t* points,
                                                  real_t* subdomain_extents,
                                                  real_t delta);

#endif
