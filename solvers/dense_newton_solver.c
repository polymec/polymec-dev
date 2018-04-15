// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/linear_algebra.h"
#include "core/norms.h"
#include "solvers/dense_newton_solver.h"
#include "nvector/nvector_serial.h"
#include "sunmatrix/sunmatrix_dense.h"
#include "sunlinsol/sunlinsol_dense.h"
#include "kinsol/kinsol.h"
#include "kinsol/kinsol_direct.h"

struct dense_newton_solver_t 
{
  // Nonlinear system information.
  int dim;
  void* context;
  int (*sys_func)(void* context, real_t* x, real_t* F);
  int (*sys_jac)(void* context, int N, real_t* x, real_t* F, 
                 real_t* work1, real_t* work2, real_t* J);
  void (*dtor)(void*);

  // Nonlinear solver.
  void* kinsol;

  // Work vectors.
  N_Vector x, x_scale, F_scale;
};

// Wrapper for system function.
static int eval_system_func(N_Vector X, N_Vector F, void* context)
{
  dense_newton_solver_t* solver = context;
  return solver->sys_func(context, NV_DATA_S(X), NV_DATA_S(F));
}

// Wrapper for system Jacobian.
static int eval_system_jac(N_Vector X, N_Vector F, SUNMatrix J,
                           void* context, N_Vector work1, N_Vector work2)
{
  dense_newton_solver_t* solver = context;
  return solver->sys_jac(context, (int)(SM_CONTENT_D(J)->N), NV_DATA_S(X), NV_DATA_S(F), 
                         NV_DATA_S(work1), NV_DATA_S(work2), SM_CONTENT_D(J)->data);
}

dense_newton_solver_t* dense_newton_solver_new(int dimension,
                                               void* context,
                                               dense_newton_system_func system_func,
                                               void (*context_dtor)(void*))
{
  return dense_newton_solver_new_with_jacobian(dimension, context, system_func,
                                               NULL, context_dtor);
}

dense_newton_solver_t* dense_newton_solver_new_with_jacobian(int dimension,
                                                             void* context,
                                                             dense_newton_system_func system_func,
                                                             dense_newton_jacobian_func jacobian_func,
                                                             void (*context_dtor)(void*))
{
  ASSERT(dimension > 0);
  ASSERT(system_func != NULL);
  dense_newton_solver_t* solver = polymec_malloc(sizeof(dense_newton_solver_t));
  solver->dim = dimension;
  solver->context = context;
  solver->dtor = context_dtor;
  solver->x = N_VNew_Serial(dimension);
  solver->x_scale = N_VNew_Serial(dimension);
  solver->F_scale = N_VNew_Serial(dimension);
  solver->kinsol = KINCreate();
  solver->sys_func = system_func;
  solver->sys_jac = jacobian_func;
  KINSetUserData(solver->kinsol, solver);
  KINInit(solver->kinsol, eval_system_func, solver->x);

  // Set up a dense linear solver.
  SUNMatrix J = SUNDenseMatrix(dimension, dimension);
  SUNLinearSolver ls = SUNDenseLinearSolver(solver->x, J);
  KINDlsSetLinearSolver(solver->kinsol, ls, J);

  // Do we have a Jacobian function?
  if (solver->sys_jac != NULL)
    KINDlsSetJacFn(solver->kinsol, eval_system_jac);

  // Use exact Newton?
  // KINSetMaxSetupCalls(solver->kinsol, 1);

  return solver;
}

void dense_newton_solver_free(dense_newton_solver_t* solver)
{
  N_VDestroy(solver->x);
  N_VDestroy(solver->x_scale);
  N_VDestroy(solver->F_scale);
  KINFree(&(solver->kinsol));
  if ((solver->context != NULL) && (solver->dtor != NULL))
    solver->dtor(solver->context);
  polymec_free(solver);
}

int dense_newton_solver_dimension(dense_newton_solver_t* solver)
{
  return solver->dim;
}

void dense_newton_solver_set_tolerances(dense_newton_solver_t* solver, real_t norm_tolerance, real_t step_tolerance)
{
  ASSERT(norm_tolerance > 0.0);
  ASSERT(step_tolerance > 0.0);
  KINSetFuncNormTol(solver->kinsol, norm_tolerance);
  KINSetScaledStepTol(solver->kinsol, step_tolerance);
}

void dense_newton_solver_set_max_iterations(dense_newton_solver_t* solver, int max_iterations)
{
  ASSERT(max_iterations > 0);
  KINSetNumMaxIters(solver->kinsol, max_iterations);
}

bool dense_newton_solver_solve(dense_newton_solver_t* solver, real_t* X, int* num_iterations)
{
  // Call the scaled solve.
  return dense_newton_solver_solve_scaled(solver, X, NULL, NULL, num_iterations);
}

bool dense_newton_solver_solve_scaled(dense_newton_solver_t* solver, real_t* X, real_t* x_scale, real_t* F_scale, int* num_iterations)
{
  // Set the x_scale and F_scale vectors.
  if (x_scale != NULL)
    memcpy(NV_DATA_S(solver->x_scale), x_scale, sizeof(real_t) * solver->dim);
  else
  {
    for (int i = 0; i < solver->dim; ++i)
      NV_Ith_S(solver->x_scale, i) = 1.0;
  }
  if (F_scale != NULL)
    memcpy(NV_DATA_S(solver->F_scale), F_scale, sizeof(real_t) * solver->dim);
  else
  {
    for (int i = 0; i < solver->dim; ++i)
      NV_Ith_S(solver->F_scale, i) = 1.0;
  }

  // Copy the values in X to the internal solution vector.
  memcpy(NV_DATA_S(solver->x), X, sizeof(real_t) * solver->dim);

  // Suspend the currently active floating point exceptions for now.
  polymec_suspend_fpe();

  // Solve.
  int status = KINSol(solver->kinsol, solver->x, KIN_LINESEARCH, 
                      solver->x_scale, solver->F_scale);

  // Reinstate the floating point exceptions.
  polymec_restore_fpe();

  if ((status == KIN_SUCCESS) || (status == KIN_INITIAL_GUESS_OK))
  {
    // Get the number of iterations it took.
    long num_iters;
    KINGetNumNonlinSolvIters(solver->kinsol, &num_iters);
    *num_iterations = (int)num_iters;

    // Copy the data back into X.
    memcpy(X, NV_DATA_S(solver->x), sizeof(real_t) * solver->dim);
    return true;
  }

  // Failed!
  return false;
}

static bool in_range(real_t x, real_t a, real_t b)
{
  return ((MIN(a, b) <= x) && (x <= MAX(a, b)));
}

// This version of Brent's method was taken from Wikipedia.
real_t brent_solve(real_t (*F)(void*, real_t), void* context, real_t x1, real_t x2, real_t tolerance, int max_iters)
{
  static const real_t delta = 1e-8;
  real_t a = x1, b = x2;
  real_t fa = F(context, a);
  real_t fb = F(context, b);
  if (fa * fb >= 0.0)
    polymec_error("brent_solve: Root is not bracketed by [x1, x2].");
  if (ABS(fa) < ABS(fb))
  {
    real_t tmp1 = b, tmp2 = fb;
    b = a;
    fb = fa;
    a = tmp1;
    fa = tmp2;
  }
  real_t c = a;
  bool mflag = true;
  real_t fs, d = 0.0, s;
  int num_iter = 0;
  do
  {
    real_t fc = F(context, c);
    if (reals_equal(fa, fc) && reals_equal(fb, fc))
    {
      s = a*fb*fc / ((fa-fb)*(fa-fc)) + 
          b*fa*fc / ((fb-fa)*(fb-fc)) + 
          c*fa*fb / ((fc-fa)*(fc-fb));
    }
    else
    {
      s = b - fb * (b-a)/(fb-fa);
    }
    if (!in_range(s, (3.0*a + b)/4.0, b) || 
        (mflag && (ABS(s-b) >= 0.5*ABS(b-c))) ||
        (!mflag && (ABS(s-b) >= 0.5*ABS(c-d))) || 
        (mflag && (ABS(b-c) < delta)) ||
        (!mflag && (ABS(c-d) < delta)))
    {
      s = 0.5 * (a + b);
      mflag = true;
    }
    else
      mflag = false;
    fs = F(context, s);
    d = c; // First assignment of d (protected by mflag).
    c = b;
    if (fa * fs < 0.0) 
    {
      b = s;
      fb = fs;
    }
    else
    {
      a = s;
      fa = fs;
    }
    if (ABS(fa) < ABS(fb))
    {
      real_t tmp1 = b, tmp2 = fb;
      b = a;
      fb = fa;
      a = tmp1;
      fa = tmp2;
    }
    ++num_iter;
  }
  while ((MAX(ABS(fb), ABS(fs)) > tolerance) && (num_iter < max_iters));
  return (ABS(fb) < ABS(fs)) ? b : s;
}

