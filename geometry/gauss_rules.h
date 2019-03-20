// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_GAUSS_RULES_H
#define POLYMEC_GAUSS_RULES_H

#include "core/polymec.h"

// Herein are functions for computing various Gaussian quadratures.

/// \addtogroup geometry geometry
///@{

/// This function fills the given array (size n each) with the quadrature
/// points and weights for the Gauss integration rule on the interval (-1, 1).
void get_gauss_points(int n, real_t* points, real_t* weights);

/// This function fills the given array (size n each) with the quadrature
/// points and weights for the Gauss-Legendre integration rule on the
/// interval (-1, 1).
void get_gauss_legendre_points(int n, real_t* points, real_t* weights);

/// This function fills the given array (size n each) with the quadrature
/// points and weights for the Gauss-Radau integration rule on the
/// interval (1, 1].
void get_gauss_radau_points(int n, real_t* points, real_t* weights);

/// This function fills the given array (size n each) with the quadrature
/// points and weights for the Gauss-Lobatto integration rule on the
/// interval [1, 1].
void get_gauss_lobatto_points(int n, real_t* points, real_t* weights);

///@}

#endif
