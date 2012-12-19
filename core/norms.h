#ifndef POLYMEC_NORMS_H
#define POLYMEC_NORMS_H

#include "core/polymec.h"

#ifdef __cplusplus
extern "C" {
#endif

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

#ifdef __cplusplus
}
#endif

#endif

