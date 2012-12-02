#include "core/brent.h"

#define BRENT_SIGN(a, b) (((b) >= 0.0) ? fabs(a) : -fabs(a))

#ifdef __cplusplus
extern "C" {
#endif

double brent_solve(brent_nl_func F, void* context, double x1, double x2, double tolerance, int max_iters, double* error)
{
  static double eps = 1e-8;
  double a = x1, b = x2;
  double fa = F(context, a), fb = F(context, b);
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
    fb = F(context, b);
    *error = fb;
  }
  return false;
}

#ifdef __cplusplus
}
#endif

