// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_NORMS_H
#define POLYMEC_NORMS_H

#include "core/polymec.h"

/// \addtogroup core core
///@{

/// Computes the L1 norm of the given vector.
static inline real_t l1_norm(real_t* vec, int dim)
{
  real_t sum = 0.0;
  for (int i = 0; i < dim; ++i)
    sum += ABS(vec[i]);
  return sum;
}

/// Computes the L2 norm of the given vector.
static inline real_t l2_norm(real_t* vec, int dim)
{
  real_t sum = 0.0;
  for (int i = 0; i < dim; ++i)
    sum += vec[i]*vec[i];
  return sqrt(sum);
}

/// Computes the L infinity norm of the given vector.
static inline real_t linf_norm(real_t* vec, int dim)
{
  real_t max = -REAL_MAX;
  for (int i = 0; i < dim; ++i)
    max = MAX(max, ABS(vec[i]));
  return max;
}

///@}

#endif

