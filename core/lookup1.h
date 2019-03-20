// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_LOOKUP1_H
#define POLYMEC_LOOKUP1_H

#include "core/polymec.h"

/// \addtogroup core core
///@{

/// \class lookup1
/// This type represents a single-parameter ("x-y") lookup table whose "x"
/// values are uniformly spaced between a minimum and maximum value.
typedef struct lookup1_t lookup1_t;

/// \enum lookup1_interp_t
/// Types of interpolation for single-parameter lookups.
typedef enum
{
  LOOKUP1_LINEAR,
  LOOKUP1_QUADRATIC
} lookup1_interp_t;

/// Constructs a new "x-y" lookup table with the given number of x values
/// distributed uniformly within [x_min, x_max], and corresponding (sorted)
/// y values.
/// \memberof lookup1
lookup1_t* lookup1_new(real_t x_min, real_t x_max, int num_values,
                       real_t* y_values, lookup1_interp_t interpolation);

/// Destroys the given lookup table.
/// \memberof lookup1
void lookup1_free(lookup1_t* table);

/// Returns the number of values in this table.
/// \memberof lookup1
int lookup1_num_values(lookup1_t* table);

/// Returns the type of interpolation used by this table.
/// \memberof lookup1
lookup1_interp_t lookup1_interpolation(lookup1_t* table);

/// Fills x_min and x_max with the lower and upper bounds of the "x" values
/// within the table.
/// \memberof lookup1
void lookup1_get_bounds(lookup1_t* table, real_t* x_min, real_t* x_max);

/// Fills y_min and y_max with the minimum and maximum y values within the
/// table.
/// \memberof lookup1
void lookup1_get_extreme_values(lookup1_t* table, real_t* y_min, real_t* y_max);

/// Returns the interpolated "y" value corresponding to the given "x" value
/// in the table, assuming that x falls within the bounds of the table.
/// \memberof lookup1
real_t lookup1_value(lookup1_t* table, real_t x);

///@}

#endif
