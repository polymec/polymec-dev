// Copyright 2012-2013 Jeffrey Johnson.
// 
// This file is part of Polymec, and is licensed under the Apache License, 
// Version 2.0 (the "License"); you may not use this file except in 
// compliance with the License. You may may find the text of the license in 
// the LICENSE file at the top-level source directory, or obtain a copy of 
// it at
// 
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef POLYMEC_NORMS_H
#define POLYMEC_NORMS_H

#include "core/polymec.h"

// Computes the L1 norm of the given vector.
static inline double l1_norm(double* vec, int dim)
{
  double sum = 0.0;
  for (int i = 0; i < dim; ++i)
    sum += fabs(vec[i]);
  return sum;
}

// Computes the L2 norm of the given vector.
static inline double l2_norm(double* vec, int dim)
{
  double sum = 0.0;
  for (int i = 0; i < dim; ++i)
    sum += vec[i]*vec[i];
  return sqrt(sum);
}

// Computes the L infinity norm of the given vector.
static inline double linf_norm(double* vec, int dim)
{
  double max = -FLT_MAX;
  for (int i = 0; i < dim; ++i)
    max = MAX(max, fabs(vec[i]));
  return max;
}

#endif

