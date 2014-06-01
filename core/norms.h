// Copyright (c) 2012-2014, Jeffrey N. Johnson
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this 
// list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice, 
// this list of conditions and the following disclaimer in the documentation 
// and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef POLYMEC_NORMS_H
#define POLYMEC_NORMS_H

#include "core/polymec.h"

// Computes the L1 norm of the given vector.
static inline real_t l1_norm(real_t* vec, int dim)
{
  real_t sum = 0.0;
  for (int i = 0; i < dim; ++i)
    sum += fabs(vec[i]);
  return sum;
}

// Computes the L2 norm of the given vector.
static inline real_t l2_norm(real_t* vec, int dim)
{
  real_t sum = 0.0;
  for (int i = 0; i < dim; ++i)
    sum += vec[i]*vec[i];
  return sqrt(sum);
}

// Computes the L infinity norm of the given vector.
static inline real_t linf_norm(real_t* vec, int dim)
{
  real_t max = -FLT_MAX;
  for (int i = 0; i < dim; ++i)
    max = MAX(max, fabs(vec[i]));
  return max;
}

#endif

