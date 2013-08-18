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

#include "core/newton.h"
#include "core/linear_algebra.h"
#include "core/norms.h"

#define BRENT_SIGN(a, b) (((b) >= 0.0) ? fabs(a) : -fabs(a))

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

double brent_solve(nonlinear_function_t F, void* context, double x1, double x2, double tolerance, int max_iters, double* error)
{
  static const double eps = 1e-8;
  double a = x1, b = x2;
  double fa, fb, deriv_dummy;
  F(context, a, &fa, &deriv_dummy);
  F(context, b, &fb, &deriv_dummy);
  if (((fa > 0.0) && (fb > 0.0)) || ((fa < 0.0) && (fb < 0.0)))
    polymec_error("brent_solve: Root is not bracketed by [x1, x2].");
  double fc = fb;
  for (int iter = 0; iter < max_iters; ++iter)
  {
    double c = 0.0, d = 0.0, e = 0.0;
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

// This estimates the Jacobian of a nonlinear system using finite differences.
static void estimate_jacobian(void* context, nonlinear_vector_function_t F, int dim, double* x, double* J)
{
  static const double eps = 1e-8;
  double f[dim], f0[dim];
  F(context, x, f0);
  for (int j = 0; j < dim; ++j)
  {
    double xj = x[j];
    double h = eps*fabs(xj);
    if (h == 0.0)
      h = eps;
    x[j] = xj + h; // Trick to reduce finite precision error.
    h = x[j] - xj;
    F(context, x, f);
    x[j] = xj;
    for (int i = 0; i < dim; ++i)
      J[dim*j+i] = (f[i] - f0[i]) / h;
  }
}

static double f_min(double* F, int dim)
{
  double sum = 0.0;
  for (int i = 0; i < dim; ++i)
    sum += F[i]*F[i];
  return 0.5*sum;
}

// This function provides line search capability for the multidimensional
// root finders. Lifted from Numerical Recipes.
static void line_search(nonlinear_system_t* system, 
                        double* x_old, double f_old, double* grad_f,
                        double* dir, double* x, double* F, double* f, 
                        double max_step, bool* check_solution)
{
  static const double eps = 1e-12;
  static const double ALF = 1e-4; // Ensures sufficient decrease in F.
  *check_solution = false;

  // Scale the direction vector for the step if necessary.
  int dim = system->dim;
  double dir_mag = 0.0;
  for (int i = 0; i < dim; ++i)
    dir_mag += dir[i]*dir[i];
  dir_mag = sqrt(dir_mag);
  if (dir_mag > max_step)
  {
    for (int i = 0; i < dim; ++i)
      dir[i] *= max_step/dir_mag;
  }

  double slope = 0.0;
  for (int i = 0; i < dim; ++i)
    slope += grad_f[i]*dir[i];
  if (slope >= 0.0)
    polymec_error("line_search: roundoff problem (slope = %g).", slope);

  double test = 0.0;
  for (int i = 0; i < dim; ++i)
  {
    double temp = fabs(dir[i])/MAX(fabs(x_old[i]), 1.0);
    if (temp > test) 
      test = temp;
  }

  double min_lambda = eps/test;
  double lambda = 1.0, lambda2 = 0.0;
  double f2 = 0.0;
  for (;;)
  {
    double temp_lambda;
    for (int i = 0; i < dim; ++i)
      x[i] = x_old[i] + lambda * dir[i];
    system->compute_F(system->context, x, F);
    *f = f_min(F, dim);
    if (lambda < min_lambda)
    {
      for (int i = 0; i < dim; ++i)
        x[i] = x_old[i];
      *check_solution = true; // Spurious convergence!
      return;
    }
    else if (*f <= f_old + ALF * lambda * slope) 
    {
      // Function f = 0.5 * F o F decreased sufficiently.
      return; 
    }
    else
    {
      if (lambda == 1.0)
      { 
        // First backtrack.
        temp_lambda = -slope/(2.0 * (*f - f_old - slope));
      }
      else
      {
        // Subsequent backtracks.
        double rhs1 = *f - f_old - lambda * slope;
        double rhs2 = f2 - f_old - lambda2 * slope;
        double a = (rhs1/(lambda*lambda) - rhs2/(lambda2*lambda2)) / (lambda - lambda2);
        double b = (-lambda2*rhs1/(lambda*lambda) + lambda*rhs2/(lambda2*lambda2))/(lambda - lambda2);
        if (a == 0.0)
          temp_lambda = -slope/(2.0*b);
        else
        {
          double disc = b*b - 3.0*a*slope;
          if (disc < 0.0)
            temp_lambda = 0.5*lambda;
          else if (b < 0.0)
            temp_lambda = (-b + sqrt(disc))/(3.0*a);
          else
            temp_lambda = -slope/(b + sqrt(disc));
        }
        if (temp_lambda > 0.5*lambda)
          temp_lambda = 0.5*lambda;
      }
    }
    lambda2 = lambda;
    f2 = *f;
    lambda = MAX(temp_lambda, 0.1*lambda);
  }
}

bool newton_solve_system(nonlinear_system_t* system, 
                         double* x, 
                         double tolerance, 
                         int max_iters,
                         int *num_iters)
{
  ASSERT(system->compute_F != NULL);

  // Some fixed parameters.
  static const double MAX_STEP = 100.0;

  int dim = system->dim;
  *num_iters = 0;

  // Check the initial guess to see if it's a root.
  double f, F[dim];
  {
    system->compute_F(system->context, x, F);
    f = f_min(F, dim);
    double test = 0.0;
    for (int i = 0; i < dim; ++i)
      test = MAX(test, fabs(F[i]));
    if (test < 0.1*tolerance)
      return true;
  }

  // Compute max_step for line searches.
  double max_step;
  {
    double sum = 0.0;
    for (int i = 0; i < dim; ++i)
      sum += x[i]*x[i];
    max_step = MAX_STEP * MAX(sqrt(sum), (double)dim);
  }

  // Iterate.
  double f_old, x_old[dim], grad_f[dim], J[dim*dim], p[dim];
  int ipiv[dim];
  for (int iter = 0; iter < max_iters; ++iter)
  {
    (*num_iters)++; // Tick.

    // Compute the Jacobian for the line search.
    if (system->compute_J != NULL)
      system->compute_J(system->context, system->compute_F, dim, x, J);
    else
      estimate_jacobian(system->context, system->compute_F, dim, x, J);

    // Now compute grad fmin = J * F.
    for (int i = 0; i < dim; ++i)
    {
      grad_f[i] = 0.0;
      for (int j = 0; j < dim; ++j)
        grad_f[i] += J[dim*i+j] * F[j];
    }

    // Store x and f, and set up the right hand side (p) for the linearized
    // system.
    for (int i = 0; i < dim; ++i)
    {
      x_old[i] = x[i];
      p[i] = -F[i];
    }
    f_old = f;

    // Solve the linear system, storing the solution in p.
    if (dim == 1)
      p[0] /= J[0];
    else if (dim == 2)
    {
      if (matrix2_det(J) == 0.0)
        polymec_error("newton_solve_system: det J == 0.");
      solve_2x2(J, p, p);
    }
    else if (dim == 3)
    {
      if (matrix3_det(J) == 0.0)
        polymec_error("newton_solve_system: det J == 0.");
      solve_3x3(J, p, p);
    }
    else
    {
      static int one = 1;
      int info;
      dgesv(&dim, &one, J, &dim, ipiv, p, &dim, &info);
      ASSERT(info == 0);
    }

    // Now do a line search.
    bool check_solution = false;
    line_search(system, x_old, f_old, grad_f, p, x, F, &f, max_step, 
                &check_solution);

    // Test for convergence on function values.
    {
      double test = 0.0;
      for (int i = 0; i < dim; ++i)
        test = MAX(test, fabs(F[i]));
      if (test < tolerance)
        return true; // Got it!
    }

    if (check_solution) // Spurious convergence?
    {
      double test = 0.0;
      double den = MAX(f, 0.5*dim);
      for (int i = 0; i < dim; ++i)
      {
        double temp = fabs(grad_f[i]) * MAX(fabs(x[i]), 1.0) / den;
        if (temp > test)
          test = temp;
      }
      return (test > tolerance);
    }

    // Test for convergence on dx.
    {
      double test = 0.0;
      for (int i = 0; i < dim; ++i)
      {
        double temp = (fabs(x[i] - x_old[i]))/MAX(fabs(x[i]), 1.0);
        if (temp > test)
          test = temp;
      }
      if (test < tolerance)
        return true; // Got it!
    }
  }
  return false;
}

bool broyden_solve_system(nonlinear_system_t* system,
                          double* x, 
                          double tolerance, 
                          int max_iters,
                          int *num_iters)
{
  // FIXME: This needs newer LAPACK joo-joo, like rank-1 QR updates.
  polymec_error("broyden_solve_system: Not implemented!");
#if 0
  static const double eps = 1e-12;
  static const double tol_min = 1e-6;
  double d[dim], f_old[dim], g[dim], p[dim],
         Q[dim*dim], R[dim*dim], tau[dim], work[dim], s[dim], t[dim], w[dim],
         x_old[dim], f[dim], fmin, fmin_old;

  bool check_solution = false;
  *num_iters = 0;

  // Check the initial guess to see if it's a root.
  {
    (*F)(x, f);
    fmin = l2_norm(f, dim);
    double test = 0.0;
    for (int i = 0; i < dim; ++i)
      test = MAX(test, fabs(f[i]));
    if (test < 0.1*tolerance)
      return true;
  }

  bool restart = true;
  for (int iter = 0; iter < max_iters; ++iter)
  {
    (*num_iters)++;
    if (restart)
    {
      // Construct a finite-difference estimate of the Jacobian and 
      // store it in R.
      estimate_jacobian(F, dim, context, x, R);

      // Compute the QR factorization. The R matrix is stored in the upper 
      // triangle of R.
      int info;
      dgeqrf(&dim, &dim, R, &dim, tau, work, &dim, &info);
      ASSERT(info == 0);
//      if (info < 0)
//        polymec_error("broyden_solve_system: Singular QR factorization in Jacobian.");

      // Form Q explicitly.
      memcpy(Q, R, dim*dim*sizeof(double));
      dorgqr(&dim, &dim, &dim, Q, &dim, tau, work, &dim, &info);
      ASSERT(info == 0);
    }
    else // Broyden update
    {
      for (int i = 0; i < dim; ++i) // s = dx
        s[i] = x[i] - x_old[i]; 
      for (int i = 0; i < dim; ++i) // t = R o s
      {
        double sum = 0.0;
        for (int j = i; j < dim; ++j)
          sum += R[dim*j+i]*s[j];
        t[i] = sum;
      }
      
      bool skip = true;
      for (int i = 0; i < dim; ++i) // w = dF - B o s
      {
        double sum = 0.0;
        for (int j = i; j < dim; ++j)
          sum += Qt[dim*i+j] * t[j];
        w[i] = f[i] - f_old[i] - sum;
        if (fabs(w[i]) >= eps * (fabs(f[i]) + fabs(f_old[i])))
          skip = false;
        else
          w[i] = 0.0;  // Hammer down noise.
      }

      if (!skip)
      {
        for (int i = 0; i < dim; ++i)
        {
          double sum = 0.0;
          for (int j = 0; j < dim; ++j)
            sum += Qt[dim*j+i]*w[j];
          t[i] = sum;
        }

        double den = 0.0;
        for (i = 0; i < dim; ++i)
          den += s[i]*s[i];
        for (int i = 0; i < dim; ++i)
          s[i] /= den;
        // FIXME: Update Q, R
        for (int i = 0; i < dim; ++i)
        {
          if (r[dim*i+i] == 0.0)
            polymec_error("broyden_solve_system: singular R.");
          D[i] = r[dim*i+i]; // Store diagonal of R in D.
        }
      }
    }

    for (int i = 0; i < dim; ++i)
    {
      double sum = 0.0;
      for (int j = i; j < dim; ++j)
        sum += Qt[dim*j+i] * f[j];
      g[i] = sum;
    }
    for (int i = dim-1; i >= 0; --i)
    {
      double sum = 0.0;
      for (int j = i; j <= i; ++j)
        sum += R[dim*i+j] * g[j];
      g[i] = sum;
    }
    for (int i = 0; i < dim; ++i)
    {
      x_old[i] = x[i];
      f_old[i] = f[i];
    }
    fmin_old = fmin;
    for (int i = 0; i < dim; ++i)
    {
      double sum = 0.0;
      for (int j = 0; j < dim; ++j)
        sum += Qt[dim*j+i] * f[j];
      p[i] = -sum;
    }

    // FIXME: Solve linear equations.
    line_search(F, dim, ...);

    // Did we do it?
    double test = 0.0;
    for (int i = 0; i < dim; ++i)
      test = MAX(test, fabs(f[i]));
    if (test < tolerance)
      return true;

    if (check_solution) // True if line search failed to find a new x.
    {
      if (restart) // Already tried reinitializing the Jacobian.
        return false;
      else
      {
        // Is the gradient of f zero?
        double test = 0.0;
        double den = MAX(fmin, 0.5*dim);
        for (int i = 0; i < dim; ++i)
        {
          double temp = fabs(g[i])*MAX(fabs(x[i]), 1.0)/den;
          if (temp > test) 
            test = temp;
        }
        if (test < tol_min) // Spurious convergence.
          return false;
        else
          restart = true;   // Try reinitializing the Jacobian.
      }
    }
    else // Successful step: use Broyden update for the next step.
    {
      restart = false;
      double test = 0.0;
      for (int i = 0; i < dim; ++i)
      {
        temp = (fabs(x[i] - x_old[i]))/MAX(fabs(x[i]), 1.0);
        if (temp > test) 
          test = temp;
      }
      if (test < eps)
        return true;
    }
  }
#endif
  return false;
}

