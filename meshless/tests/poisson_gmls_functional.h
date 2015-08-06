// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POISSON_GMLS_FUNCTIONAL_H
#define POISSON_GMLS_FUNCTIONAL_H

#include "core/point_cloud.h"
#include "meshless/gmls_functional.h"

// Constructs a GMLS functional for solving Poisson's equation on the 
// interior of a domain, using a polynomial basis of the given degree.
gmls_functional_t* poisson_gmls_functional_new(int degree,
                                               point_cloud_t* points,
                                               real_t* subdomain_extents);

// Constructs a GMLS functional for solving Poisson's equation on a Dirichlet 
// boundary, using a polynomial basis of the given degree.
gmls_functional_t* poisson_dirichlet_gmls_functional_new(int degree,
                                                         point_cloud_t* points,
                                                         real_t* subdomain_extents);

#endif
