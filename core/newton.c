#include "core/newton.h"

#ifdef __cplusplus
extern "C" {
#endif

bool newton_solve(nonlinear_function_t F, void* context, double* x, double min, double max, double tolerance, int max_iters)
{
  double f, dfdx;
  F(context, *x, &f, &dfdx);
  if (f < tolerance) 
    return true;
  for (int iter = 0; iter < max_iters; ++iter)
  {
    *x -= f/dfdx;
    F(context, *x, &f, &dfdx);
    if (f < tolerance) 
      return true;
  }
  return false;
}

bool picard_solve(nonlinear_function_t F, void* context, double* x, double min, double max, double tolerance, int max_iters)
{
  double f, dfdx;
  F(context, *x, &f, &dfdx);
  if (f < tolerance) 
    return true;
  for (int iter = 0; iter < max_iters; ++iter)
  {
    *x = f;
    F(context, *x, &f, &dfdx);
    if (f < tolerance) 
      return true;
  }
  return false;
}

bool newton_picard_solve(nonlinear_function_t F, void* context, double* x, double min, double max, double tolerance, int newton_iters, int picard_iters, int num_cycles)
{
  for (int cycle = 0; cycle < num_cycles; ++cycle)
  {
    if (newton_solve(F, context, x, min, max, newton_iters, tolerance))
      return true;
    if (picard_solve(F, context, x, min, max, picard_iters, tolerance))
      return true;
  }
  return false;
}

#ifdef __cplusplus
}
#endif

