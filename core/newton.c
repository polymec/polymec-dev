// Copyright (c) 2012-2013, Jeffrey N. Johnson
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

#include "core/newton.h"
#include "core/linear_algebra.h"
#include "core/norms.h"
#include "kinsol/kinsol.h"
#include "kinsol/kinsol_dense.h"

struct newton_solver_t 
{
  // Nonlinear system information.
  int dim;
  void* context;
  void (*dtor)(void*);

  // Nonlinear solver.
  void* kinsol;

  // Work vectors.
  N_Vector x, x_scale, F_scale;
};

newton_solver_t* newton_solver_new(int dimension,
                                   void* context,
                                   newton_system_func_t system_func,
                                   void (*context_dtor)(void*))
{
  ASSERT(dimension > 0);
  ASSERT(system_func != NULL);
  newton_solver_t* solver = malloc(sizeof(newton_solver_t));
  solver->dim = dimension;
  solver->context = context;
  solver->dtor = context_dtor;
  solver->x = N_VNew_Serial(dimension);
  solver->x_scale = N_VNew_Serial(dimension);
  solver->F_scale = N_VNew_Serial(dimension);
  solver->kinsol = KINCreate();
  KINInit(solver->kinsol, system_func, solver->x);
  KINDense(solver->kinsol, dimension);

  // Use exact Newton?
  // KINSetMaxSetupCalls(solver->kinsol, 1);

  return solver;
}

void newton_solver_free(newton_solver_t* solver)
{
  N_VDestroy(solver->x);
  N_VDestroy(solver->x_scale);
  N_VDestroy(solver->F_scale);
  KINFree(solver->kinsol);
  if ((solver->context != NULL) && (solver->dtor != NULL))
    solver->dtor(solver->context);
  free(solver);
}

int newton_solver_dimension(newton_solver_t* solver)
{
  return solver->dim;
}

void newton_solver_set_tolerances(newton_solver_t* solver, double norm_tolerance, double step_tolerance)
{
  ASSERT(norm_tolerance > 0.0);
  ASSERT(step_tolerance > 0.0);
  KINSetFuncNormTol(solver->kinsol, norm_tolerance);
  KINSetScaledStepTol(solver->kinsol, step_tolerance);
}

void newton_solver_set_max_iterations(newton_solver_t* solver, int max_iterations)
{
  ASSERT(max_iterations > 0);
  KINSetNumMaxIters(solver->kinsol, max_iterations);
}

bool newton_solver_solve(newton_solver_t* solver, double* X, int* num_iterations)
{
  // Call the scaled solve.
  return newton_solver_solve_scaled(solver, X, NULL, NULL, num_iterations);
}

bool newton_solver_solve_scaled(newton_solver_t* solver, double* X, double* x_scale, double* F_scale, int* num_iterations)
{
  // Set the x_scale and F_scale vectors.
  if (x_scale != NULL)
    memcpy(NV_DATA_S(solver->x_scale), x_scale, sizeof(double) * solver->dim);
  else
  {
    for (int i = 0; i < solver->dim; ++i)
      NV_Ith_S(solver->x_scale, i) = 1.0;
  }
  if (F_scale != NULL)
    memcpy(NV_DATA_S(solver->F_scale), F_scale, sizeof(double) * solver->dim);
  else
  {
    for (int i = 0; i < solver->dim; ++i)
      NV_Ith_S(solver->F_scale, i) = 1.0;
  }

  // Copy the values in X to the internal solution vector.
  memcpy(NV_DATA_S(solver->x), X, sizeof(double) * solver->dim);

  // Solve.
  int status = KINSol(solver->kinsol, solver->x, KIN_LINESEARCH, 
                      solver->x_scale, solver->F_scale);
  if ((status == KIN_SUCCESS) || (status == KIN_INITIAL_GUESS_OK))
  {
    // Get the number of iterations it took.
    long num_iters;
    KINGetNumNonlinSolvIters(solver->kinsol, &num_iters);
    *num_iterations = (int)num_iters;

    // Copy the data back into X.
    memcpy(X, NV_DATA_S(solver->x), sizeof(double) * solver->dim);
    return true;
  }

  // Failed!
  return false;
}

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

