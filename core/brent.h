#ifndef POLYMEC_BRENT_H
#define POLYMEC_BRENT_H

#include "core/polymec.h"

#ifdef __cplusplus
extern "C" {
#endif

// This defines a function that computes the value of a 
// (real-valued) function of a single variable given a context.
typedef double (*brent_nl_func)(void*, double);

// Solve the nonlinear equation F(x) = 0 using Brent's method on 
// the interval [x1, x2] with the specified maximum number of iterations 
// and desired tolerance. Returns the solution x.
double brent_solve(brent_nl_func F, void* context, double x1, double x2, double tolerance, int max_iters, double* error); 

#ifdef __cplusplus
}
#endif

#endif

