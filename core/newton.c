#include "core/newton.h"

#define BRENT_SIGN(a, b) (((b) >= 0.0) ? fabs(a) : -fabs(a))

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

bool picard_solve(nonlinear_function_t F, void* context, double* x, double tolerance, int max_iters)
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

bool picard_solve_system(int dim, nonlinear_vector_function_t F, void* context, double* x, double tolerance, int max_iters)
{
  double f[dim], J[dim*dim];
  F(context, *x, f, J);
  double l2f = 0.0;
  for (int i = 0; i < dim; ++i)
    l2f += (f[i]*f[i]);
  if (l2f < tolerance) 
    return true;
  for (int iter = 0; iter < max_iters; ++iter)
  {
    for (int i = 0; i < dim; ++i)
      x[i] = f[i];
    F(context, *x, f, J);
    l2f = 0.0;
    for (int i = 0; i < dim; ++i)
      l2f += (f[i]*f[i]);
    if (l2f < tolerance) 
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
    if (picard_solve(F, context, x, picard_iters, tolerance))
      return true;
  }
  return false;
}

double brent_solve(nonlinear_function_t F, void* context, double x1, double x2, double tolerance, int max_iters, double* error)
{
  static double eps = 1e-8;
  double a = x1, b = x2;
  double fa, fb, deriv_dummy;
  F(context, a, &fb, &deriv_dummy);
  F(context, b, &fb, &deriv_dummy);
  if (((fa > 0.0) && (fb > 0.0)) || ((fa < 0.0) && (fb < 0.0)))
    polymec_error("brent_solve: Root is not bracketed by [x1, x2].");
  double fc = fb;
  for (int iter = 0; iter < max_iters; ++iter)
  {
    double c, d, e;
    if (((fb > 0.0) && (fc > 0.0)) || (((fb < 0.0) && (fc < 0.0))))
    {
      c = a;
      fc = fa;
      e = d = b - a;
    }
    if (fabs(fc) < fabs(fb)) 
    {
      c = a;
      b = c;
      c = a;
      fa = fb;
      fb = fc;
      fc = fa;
    }
    double tol1 = 2.0*eps*fabs(b)+0.5*tolerance;
    double xm = 0.5*(c - b);
    if ((fabs(xm) <= tol1) || (fb == 0.0))
      return b;
    if ((fabs(e) >= tol1) && (fabs(fa) > fabs(fb))) 
    {
      double p, q, r, s;
      s = fb/fa; // Attempt inverse quadratic interpolation
      if (a == c)
      {
        p = 2.0*xm*s;
        q = 1.0 - s;
      }
      else
      {
        q = fa/fc;
        r = fb/fc;
        p = s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
        q = (q-1.0)*(r-1.0)*(s-1.0);
      }
      if (p > 0.0) 
        q = -q;
      p = fabs(p);
      double min1 = 3.0*xm*q - fabs(tol1*q);
      double min2 = fabs(e*q);
      if (2.0*p < ((min1 < min2) ? min1 : min2))
      {
        e = d;
        d = p/q;
      }
      else
      {
        d = xm;
        e = d;
      }
    }
    else
    {
      d = xm;
      e = d;
    }
    a = b;
    fa = fb;
    if (fabs(d) > tol1)
      b += d;
    else
      b += BRENT_SIGN(tol1, xm);
    F(context, b, &fb, &deriv_dummy);
    *error = fb;
  }
  return false;
}

#ifdef __cplusplus
}
#endif

