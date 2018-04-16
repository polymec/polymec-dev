// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/sundials_helpers.h"
#include "core/timer.h"
#include "solvers/dae_solver.h"

#include "ida/ida.h"
#include "ida/ida_spils.h"

// JFNK stuff.
#include "sunlinsol/sunlinsol_spgmr.h"
#include "sunlinsol/sunlinsol_spfgmr.h"
#include "sunlinsol/sunlinsol_spbcgs.h"
#include "sunlinsol/sunlinsol_sptfqmr.h"

// Stuff for generalized solvers.
#include "ida/ida_impl.h"

// These symbols represent special short-hand arguments for equation types 
// and constraints.
static dae_equation_t DAE_ALL_ALGEBRAIC_LOC = DAE_ALGEBRAIC;
dae_equation_t* DAE_ALL_ALGEBRAIC = &DAE_ALL_ALGEBRAIC_LOC;

static dae_equation_t DAE_ALL_DIFFERENTIAL_LOC = DAE_DIFFERENTIAL;
dae_equation_t* DAE_ALL_DIFFERENTIAL = &DAE_ALL_DIFFERENTIAL_LOC;

static dae_constraint_t DAE_ALL_UNCONSTRAINED_LOC = DAE_UNCONSTRAINED;
dae_constraint_t* DAE_ALL_UNCONSTRAINED = &DAE_ALL_UNCONSTRAINED_LOC;

static dae_constraint_t DAE_ALL_NEGATIVE_LOC = DAE_NEGATIVE;
dae_constraint_t* DAE_ALL_NEGATIVE = &DAE_ALL_NEGATIVE_LOC;

static dae_constraint_t DAE_ALL_NONPOSITIVE_LOC = DAE_NONPOSITIVE;
dae_constraint_t* DAE_ALL_NONPOSITIVE = &DAE_ALL_NONPOSITIVE_LOC;

static dae_constraint_t DAE_ALL_NONNEGATIVE_LOC = DAE_NONNEGATIVE;
dae_constraint_t* DAE_ALL_NONNEGATIVE = &DAE_ALL_NONNEGATIVE_LOC;

static dae_constraint_t DAE_ALL_POSITIVE_LOC = DAE_POSITIVE;
dae_constraint_t* DAE_ALL_POSITIVE = &DAE_ALL_POSITIVE_LOC;

struct dae_solver_t 
{
  char* name;
  MPI_Comm comm;
  int num_local_values, num_remote_values;
  void* context;
  int order;
  bool initialized;
  int (*F)(void* context, real_t t, real_t* U, real_t* U_dot, real_t* F);
  void (*dtor)(void* context);

  // IDA data structures.
  void* ida;
  N_Vector U, U_dot;
  real_t* U_with_ghosts;
  real_t* U_dot_with_ghosts;
  real_t t;
  char* status_message; // status of most recent integration.
  real_t max_dt, stop_time;

  // JFNK stuff.
  int max_krylov_dim;
  newton_pc_t* precond;
  jfnk_dae_krylov_t solver_type;
  int (*Jy)(void* context, real_t t, real_t* U, real_t alpha, real_t* U_dot, real_t* F,
            real_t* y, real_t* Jy, real_t* tmp1, real_t* tmp2);
  SUNLinearSolver ls;

  // Generalized adaptor stuff.
  real_t sqrtN;
  int (*reset_func)(void* context, real_t t, real_t* U, real_t* U_dot);
  int (*setup_func)(void* context, 
                    real_t alpha, 
                    int step,
                    real_t t, 
                    real_t* U_pred, 
                    real_t* U_dot_pred, 
                    real_t* F_pred, 
                    real_t* work1, real_t* work2, real_t* work3);
  int (*solve_func)(void* context, 
                    real_t t, 
                    real_t* U,
                    real_t* U_dot,
                    real_t* F,
                    real_t* W, 
                    real_t res_norm_tol,
                    real_t* B,
                    int* num_iters);
  int num_linear_iterations, num_linear_conv_failures;

  // Error weight function.
  void (*compute_weights)(void* context, real_t* y, real_t* weights);
  real_t* error_weights;
};

static char* get_status_message(int status, real_t current_time)
{
  char* status_message = NULL;
  if (status == IDA_TOO_MUCH_WORK)
  {
    char err[1024];
    snprintf(err, 1024, "solver stopped at t = %g after maximum number of steps.", current_time);
    status_message = string_dup(err);
  }
  else if (status == IDA_TOO_MUCH_ACC)
    status_message = string_dup("solver could not achieve desired level of accuracy.");
  else if (status == IDA_ERR_FAIL)
    status_message = string_dup("solver encountered too many error test failures.");
  else if (status == IDA_CONV_FAIL)
    status_message = string_dup("solver encountered too many convergence test failures.");
  else if (status == IDA_LINIT_FAIL)
    status_message = string_dup("solver's linear solver failed to initialize.");
  else if (status == IDA_LSETUP_FAIL)
    status_message = string_dup("solver's linear solver setup failed.");
  else if (status == IDA_LSOLVE_FAIL)
    status_message = string_dup("solver's linear solver failed.");
  else if (status == IDA_CONSTR_FAIL)
    status_message = string_dup("solver's inequality constraints were violated unrecoverably.");
  else if (status == IDA_RES_FAIL)
    status_message = string_dup("solver's residual function failed unrecoverably.");
  else if (status == IDA_REP_RES_ERR)
    status_message = string_dup("solver encountered too many recoverable RHS failures.");
  else if (status == IDA_RTFUNC_FAIL)
    status_message = string_dup("solver encountered a failure in the rootfinding function.");
  return status_message;
}

// This function wraps around the user-supplied right hand side.
static int evaluate_residual(real_t t, N_Vector U, N_Vector U_dot, 
                             N_Vector F, void* context)
{
  dae_solver_t* integ = context;
  real_t* xx = NV_DATA(U);
  real_t* xxd = NV_DATA(U_dot);
  real_t* Fx = NV_DATA(F);

  // Evaluate the residual using vectors with ghosts.
  memcpy(integ->U_with_ghosts, xx, sizeof(real_t) * integ->num_local_values);
  memcpy(integ->U_dot_with_ghosts, xxd, sizeof(real_t) * integ->num_local_values);
  return integ->F(integ->context, t, integ->U_with_ghosts, integ->U_dot_with_ghosts, Fx);
}

// This function sets up the preconditioner data within the solver.
static int set_up_preconditioner(real_t t, N_Vector U, N_Vector U_dot, N_Vector F,
                                 real_t cj, void* context)
{
  START_FUNCTION_TIMER();
  dae_solver_t* integ = context;

  // Form the linear combination dFdx + cj * dFd(U_dot).
  memcpy(integ->U_with_ghosts, NV_DATA(U), sizeof(real_t) * integ->num_local_values);
  memcpy(integ->U_dot_with_ghosts, NV_DATA(U_dot), sizeof(real_t) * integ->num_local_values);
  newton_pc_setup(integ->precond, 0.0, 1.0, cj, t, integ->U_with_ghosts, integ->U_dot_with_ghosts);
  STOP_FUNCTION_TIMER();
  return 0;
}

// This function solves the preconditioner equation. On input, the vector r 
// contains the right-hand side of the preconditioner system, and on output 
// it contains the solution to the system.
static int solve_preconditioner_system(real_t t, N_Vector U, N_Vector U_dot, 
                                       N_Vector F, N_Vector r, N_Vector z, 
                                       real_t cj, real_t delta, void* context)
{
  START_FUNCTION_TIMER();
  dae_solver_t* integ = context;
  
  // FIXME: Apply scaling if needed.

  // Solve it.
  int result = (newton_pc_solve(integ->precond, t, NV_DATA(U), NV_DATA(U_dot),
                                NV_DATA(r), NV_DATA(z))) ? 0 : 1;
  STOP_FUNCTION_TIMER();
  return result;
}

// Adaptor for J*y setup function
static int set_up_Jy(real_t t, N_Vector y, N_Vector ydot, N_Vector Jy, real_t cj, void* context)
{
  // Nothing here!
  return 0;
}

// Adaptor for J*y function.
static int eval_Jy(real_t tt, N_Vector uu, N_Vector up, N_Vector ff,
                   N_Vector y, N_Vector Jy, real_t cj, 
                   void *context, N_Vector tmp1, N_Vector tmp2)
{
  START_FUNCTION_TIMER();
  dae_solver_t* integ = context;
  real_t* U = NV_DATA(uu);
  real_t* Udot = NV_DATA(up);
  real_t* F = NV_DATA(ff);
  real_t* yy = NV_DATA(y);
  real_t* Jyy = NV_DATA(Jy);
  real_t* tmp11 = NV_DATA(tmp1);
  real_t* tmp22 = NV_DATA(tmp2);

  // Make sure we use ghosts.
  memcpy(integ->U_with_ghosts, U, sizeof(real_t) * integ->num_local_values);
  memcpy(integ->U_dot_with_ghosts, Udot, sizeof(real_t) * integ->num_local_values);
  int status = integ->Jy(integ->context, tt, U, cj, Udot, F, yy, Jyy, tmp11, tmp22);
  STOP_FUNCTION_TIMER();
  return status;
}

// We use floating point comparisons to assign inequality relations in this 
// function, so we disable related compiler warnings.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wfloat-equal"
static void set_up_equations_and_constraints(dae_solver_t* integ, 
                                             dae_equation_t* equation_types,
                                             dae_constraint_t* constraints)
{
  // We set the equation types.
  N_Vector array = N_VNew(integ->comm, integ->num_local_values);
  bool has_algebraic_eqns = false;
  if (equation_types == DAE_ALL_ALGEBRAIC)
  {
    has_algebraic_eqns = true;
    for (int i = 0; i < integ->num_local_values; ++i)
      NV_Ith(array, i) = 0.0;
  }
  else if (equation_types == DAE_ALL_DIFFERENTIAL)
  {
    for (int i = 0; i < integ->num_local_values; ++i)
      NV_Ith(array, i) = 1.0;
  }
  else
  {
    for (int i = 0; i < integ->num_local_values; ++i)
    {
      if (equation_types[i] == DAE_ALGEBRAIC)
      {
        has_algebraic_eqns = true;
        NV_Ith(array, i) = 0.0;
      }
      else
        NV_Ith(array, i) = 1.0;
    }
  }
  IDASetId(integ->ida, array);

  // Set up constraints.
  if (constraints != DAE_ALL_UNCONSTRAINED)
  {
    if (constraints == DAE_ALL_NEGATIVE)
    {
      for (int i = 0; i < integ->num_local_values; ++i)
        NV_Ith(array, i) = -2.0;
    }
    else if (constraints == DAE_ALL_NONPOSITIVE)
    {
      for (int i = 0; i < integ->num_local_values; ++i)
        NV_Ith(array, i) = -1.0;
    }
    else if (constraints == DAE_ALL_NONNEGATIVE)
    {
      for (int i = 0; i < integ->num_local_values; ++i)
        NV_Ith(array, i) = 1.0;
    }
    else if (constraints == DAE_ALL_POSITIVE)
    {
      for (int i = 0; i < integ->num_local_values; ++i)
        NV_Ith(array, i) = 2.0;
    }
    else
    {
      for (int i = 0; i < integ->num_local_values; ++i)
      {
        real_t constraint = constraints[i];
        if (constraint == DAE_UNCONSTRAINED)
          NV_Ith(array, i) = 0.0;
        else if (constraint == DAE_NEGATIVE)
          NV_Ith(array, i) = -2.0;
        else if (constraint == DAE_NONPOSITIVE)
          NV_Ith(array, i) = -1.0;
        else if (constraint == DAE_NONNEGATIVE)
          NV_Ith(array, i) = 1.0;
        else // (constraint == DAE_POSITIVE)
          NV_Ith(array, i) = 2.0;
      }
    }
    IDASetConstraints(integ->ida, array);
  }
  N_VDestroy(array);

  // If we have algebraic equations, we suppress their contributions to 
  // the estimate of the truncation error.
  if (has_algebraic_eqns)
    IDASetSuppressAlg(integ->ida, 1);
}
#pragma GCC diagnostic pop
#pragma clang diagnostic pop

dae_solver_t* jfnk_dae_solver_new(int order,
                                  MPI_Comm comm,
                                  dae_equation_t* equation_types,
                                  dae_constraint_t* constraints,
                                  int num_local_values,
                                  int num_remote_values,
                                  void* context,
                                  int (*F_func)(void* context, real_t t, real_t* U, real_t* U_dot, real_t* F),
                                  int (*Jy_func)(void* context, real_t t, real_t* U, real_t alpha, real_t* U_dot, real_t* F,
                                                 real_t* y, real_t* Jy, real_t* tmp1, real_t* tmp2),
                                  void (*dtor)(void* context),
                                  newton_pc_t* precond,
                                  jfnk_dae_krylov_t solver_type,
                                  int max_krylov_dim)
{
  ASSERT(order > 0);
  ASSERT(order <= 5);
  ASSERT(equation_types != NULL);
  ASSERT(constraints != NULL);
  ASSERT(num_local_values > 0);
  ASSERT(num_remote_values >= 0);
  ASSERT(F_func != NULL);
  ASSERT(precond != NULL);
  ASSERT(newton_pc_side(precond) == NEWTON_PC_LEFT); // left PC only!
  ASSERT(max_krylov_dim >= 3);

  dae_solver_t* integ = polymec_malloc(sizeof(dae_solver_t));
  integ->name = string_dup("JFNK Differential-Algebraic solver");
  integ->context = context;
  integ->comm = comm;
  integ->order = order;
  integ->solver_type = solver_type;
  integ->t = 0.0;
  integ->F = F_func;
  integ->Jy = Jy_func;
  integ->dtor = dtor;
  integ->num_local_values = num_local_values;
  integ->num_remote_values = num_remote_values;
  integ->precond = precond;
  integ->max_krylov_dim = max_krylov_dim;
  integ->initialized = false;
  integ->max_dt = REAL_MAX;
  integ->status_message = NULL;
  integ->error_weights = NULL;

  integ->reset_func = NULL;
  integ->setup_func = NULL;
  integ->solve_func = NULL;

  // Set up IDA and accessories.
  integ->U = N_VNew(comm, num_local_values);
  integ->U_with_ghosts = polymec_malloc(sizeof(real_t) * (num_local_values + num_remote_values));
  integ->U_dot = N_VNew(comm, num_local_values);
  integ->U_dot_with_ghosts = polymec_malloc(sizeof(real_t) * (num_local_values + num_remote_values));
  memset(NV_DATA(integ->U), 0, sizeof(real_t) * num_local_values);
  memset(NV_DATA(integ->U_dot), 0, sizeof(real_t) * num_local_values);

  integ->ida = IDACreate();
  IDASetMaxOrd(integ->ida, integ->order);
  IDASetUserData(integ->ida, integ);
  IDAInit(integ->ida, evaluate_residual, integ->t, integ->U, integ->U_dot);

  // Select the particular type of Krylov method for the underlying linear solves.
  if (solver_type == JFNK_DAE_GMRES)
  {
    integ->ls = SUNSPGMR(integ->U, PREC_LEFT, max_krylov_dim);
    // We use modified Gram-Schmidt orthogonalization.
    SUNSPGMRSetGSType(integ->ls, MODIFIED_GS);
    IDASpilsSetLinearSolver(integ->ida, integ->ls);
  }
  else if (solver_type == JFNK_DAE_BICGSTAB)
  {
    integ->ls = SUNSPBCGS(integ->U, PREC_LEFT, max_krylov_dim);
    IDASpilsSetLinearSolver(integ->ida, integ->ls);
  }
  else
  {
    integ->ls = SUNSPTFQMR(integ->U, PREC_LEFT, max_krylov_dim);
    IDASpilsSetLinearSolver(integ->ida, integ->ls);
  }

  // Set up the equation types and constraints.
  set_up_equations_and_constraints(integ, equation_types, constraints);

  // Set up the Jacobian function if given.
  if (integ->Jy != NULL)
    IDASpilsSetJacTimes(integ->ida, set_up_Jy, eval_Jy);

  // Set up preconditioner machinery.
  IDASpilsSetPreconditioner(integ->ida, set_up_preconditioner,
                            solve_preconditioner_system);

  // Set some default tolerances:
  // relative error of 1e-4 means errors are controlled to 0.01%.
  // absolute error is set to 1 because it's completely problem dependent.
  dae_solver_set_tolerances(integ, 1e-4, 1.0);

  // Set up a maximum number of steps to take during the integration.
//  IDASetMaxNumSteps(integ->ida, 500); // default is 500.

  return integ;
}

static int dae_linit(IDAMem ida_mem)
{
  dae_solver_t* dae = ida_mem->ida_user_data;
  real_t t = ida_mem->ida_tn;
  real_t* U = NV_DATA(ida_mem->ida_yy);
  real_t* U_dot = NV_DATA(ida_mem->ida_yp);
  dae->sqrtN = -1.0;
  dae->num_linear_iterations = 0;
  dae->num_linear_conv_failures = 0;
  return dae->reset_func(dae->context, t, U, U_dot);
}

static int dae_lsetup(IDAMem ida_mem, N_Vector yyp, N_Vector ypp, N_Vector resp,
                      N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3)
{
  dae_solver_t* dae = ida_mem->ida_user_data;
  real_t t = ida_mem->ida_tn;
  real_t alpha = ida_mem->ida_cj;
  int step = (int)ida_mem->ida_nst;
  real_t* U_pred = NV_DATA(yyp);
  real_t* U_dot_pred = NV_DATA(ypp);
  real_t* F_pred = NV_DATA(resp);
  real_t* work1 = NV_DATA(vtemp1);
  real_t* work2 = NV_DATA(vtemp2);
  real_t* work3 = NV_DATA(vtemp3);
  if (dae->sqrtN <= 0.0)
  {
    N_VConst(1.0, vtemp1);
    dae->sqrtN = sqrt(N_VDotProd(vtemp1, vtemp1));
  }
  return dae->setup_func(dae->context, alpha, step, t, U_pred, U_dot_pred, F_pred,
                         work1, work2, work3);
}

static int dae_lsolve(IDAMem ida_mem, N_Vector b, N_Vector weight, N_Vector ycur,
                      N_Vector ypcur, N_Vector rescur)
{
  dae_solver_t* dae = ida_mem->ida_user_data;
  real_t t = ida_mem->ida_tn;
  real_t* U = NV_DATA(ycur);
  real_t* U_dot = NV_DATA(ypcur);
  real_t* F = NV_DATA(rescur);
  real_t* W = NV_DATA(weight);
  real_t res_norm_tol = 0.5 * dae->sqrtN * ida_mem->ida_epsNewt;
  real_t* B = NV_DATA(b);
  int num_iters = 0;
  int result = dae->solve_func(dae->context, t, U, U_dot, F, W, res_norm_tol, B, &num_iters);
  dae->num_linear_iterations += num_iters;
  if (result != 0)
    ++(dae->num_linear_conv_failures);
  return result;
}

static int dae_lperf(IDAMem ida_mem, int perftask)
{
  return 0;
}

static int dae_lfree(IDAMem ida_mem)
{
  return 0;
}

dae_solver_t* dae_solver_new(const char* name,
                             int order, 
                             MPI_Comm comm,
                             dae_equation_t* equation_types,
                             dae_constraint_t* constraints,
                             int num_local_values, 
                             int num_remote_values, 
                             void* context, 
                             int (*F_func)(void* context, real_t t, real_t* U, real_t* U_dot, real_t* F),
                             int (*reset_func)(void* context, real_t t, real_t* U, real_t* U_dot),
                             int (*setup_func)(void* context, 
                                               real_t alpha, 
                                               int step,
                                               real_t t, 
                                               real_t* U_pred, 
                                               real_t* U_dot_pred, 
                                               real_t* F_pred, 
                                               real_t* work1, real_t* work2, real_t* work3),
                             int (*solve_func)(void* context, 
                                               real_t t, 
                                               real_t* U,
                                               real_t* U_dot,
                                               real_t* F,
                                               real_t* W, 
                                               real_t res_norm_tol,
                                               real_t* B,
                                               int* num_iters), 
                             void (*dtor)(void* context))
{
  ASSERT(order > 0);
  ASSERT(order <= 5);
  ASSERT(equation_types != NULL);
  ASSERT(constraints != NULL);
  ASSERT(num_local_values > 0);
  ASSERT(num_remote_values >= 0);
  ASSERT(F_func != NULL);

  dae_solver_t* integ = polymec_malloc(sizeof(dae_solver_t));
  integ->name = string_dup(name);
  integ->context = context;
  integ->comm = comm;
  integ->order = order;
  integ->t = 0.0;
  integ->F = F_func;
  integ->Jy = NULL;
  integ->dtor = dtor;
  integ->num_local_values = num_local_values;
  integ->num_remote_values = num_remote_values;
  integ->precond = NULL;
  integ->max_krylov_dim = -1;
  integ->initialized = false;
  integ->max_dt = REAL_MAX;
  integ->status_message = NULL;
  integ->error_weights = NULL;

  // Set up IDA and accessories.
  integ->U = N_VNew(comm, num_local_values);
  integ->U_with_ghosts = polymec_malloc(sizeof(real_t) * (num_local_values + num_remote_values));
  integ->U_dot = N_VNew(comm, num_local_values);
  integ->U_dot_with_ghosts = polymec_malloc(sizeof(real_t) * (num_local_values + num_remote_values));
  memset(NV_DATA(integ->U), 0, sizeof(real_t) * num_local_values);
  memset(NV_DATA(integ->U_dot), 0, sizeof(real_t) * num_local_values);

  integ->ida = IDACreate();
  IDASetMaxOrd(integ->ida, integ->order);
  IDASetUserData(integ->ida, integ);
  IDAInit(integ->ida, evaluate_residual, integ->t, integ->U, integ->U_dot);

  // Set up the INK solver.
  integ->reset_func = reset_func;
  integ->setup_func = setup_func;
  integ->solve_func = solve_func;
  IDAMem ida_mem = integ->ida;
  ida_mem->ida_linit = dae_linit;
  ida_mem->ida_lsetup = dae_lsetup;
  ida_mem->ida_lsolve = dae_lsolve;
  ida_mem->ida_lperf = dae_lperf;
  ida_mem->ida_lfree = dae_lfree;

  // Set up the equation types and constraints.
  set_up_equations_and_constraints(integ, equation_types, constraints);

  // Set some default tolerances:
  // relative error of 1e-4 means errors are controlled to 0.01%.
  // absolute error is set to 1 because it's completely problem dependent.
  dae_solver_set_tolerances(integ, 1e-4, 1.0);

  return integ;
}

void dae_solver_free(dae_solver_t* integ)
{
  string_free(integ->name);

  // Kill the preconditioner stuff.
  if (integ->precond != NULL)
    newton_pc_free(integ->precond);

  // Kill the IDA stuff.
  N_VDestroy(integ->U_dot);
  polymec_free(integ->U_with_ghosts);
  N_VDestroy(integ->U);
  polymec_free(integ->U_dot_with_ghosts);
  SUNLinSolFree(integ->ls);
  IDAFree(&integ->ida);

  // Kill the rest.
  if (integ->status_message != NULL)
    polymec_free(integ->status_message);
  if ((integ->context != NULL) && (integ->dtor != NULL))
    integ->dtor(integ->context);
  if (integ->error_weights != NULL)
    polymec_free(integ->error_weights);
  polymec_free(integ);
}

char* dae_solver_name(dae_solver_t* integ)
{
  return integ->name;
}

void* dae_solver_context(dae_solver_t* integ)
{
  return integ->context;
}

int dae_solver_order(dae_solver_t* integ)
{
  return integ->order;
}

newton_pc_t* dae_solver_preconditioner(dae_solver_t* solver)
{
  return solver->precond;
}

void dae_solver_set_tolerances(dae_solver_t* integ,
                               real_t relative_tol, real_t absolute_tol)
{
  ASSERT(relative_tol > 0.0);
  ASSERT(absolute_tol > 0.0);

  // Clear any existing error weight function.
  integ->compute_weights = NULL;
  if (integ->error_weights != NULL)
  {
    polymec_free(integ->error_weights);
    integ->error_weights = NULL;
  }


  // Set the tolerances.
  IDASStolerances(integ->ida, relative_tol, absolute_tol);
}

// Constant error weight adaptor function.
static void use_constant_weights(void* context, real_t* y, real_t* weights)
{
  dae_solver_t* integ = context;
  ASSERT(integ->error_weights != NULL);
  memcpy(weights, integ->error_weights, sizeof(real_t) * integ->num_local_values);
}

void dae_solver_set_error_weights(dae_solver_t* integ, real_t* weights)
{
#ifndef NDEBUG
  // Check for non-negativity and total positivity.
  real_t total = 0.0;
  for (int i = 0; i < integ->num_local_values; ++i)
  {
    ASSERT(weights[i] >= 0.0);
    total += weights[i];
  }
  ASSERT(total > 0.0);
#endif

  if (integ->error_weights == NULL)
    integ->error_weights = polymec_malloc(sizeof(real_t) * integ->num_local_values);
  memcpy(integ->error_weights, weights, sizeof(real_t) * integ->num_local_values);
  dae_solver_set_error_weight_function(integ, use_constant_weights);
}

// Error weight adaptor function.
static int compute_error_weights(N_Vector y, N_Vector ewt, void* context)
{
  dae_solver_t* integ = context;
  integ->compute_weights(integ->context, NV_DATA(y), NV_DATA(ewt));
  return 0;
}

void dae_solver_set_error_weight_function(dae_solver_t* integ,
                                          void (*compute_weights)(void* context, real_t* y, real_t* weights))
{
  ASSERT(compute_weights != NULL);
  integ->compute_weights = compute_weights;
  IDAWFtolerances(integ->ida, compute_error_weights);
}

void dae_solver_eval_residual(dae_solver_t* integ, real_t t, real_t* U, real_t* U_dot, real_t* F);
void dae_solver_eval_residual(dae_solver_t* integ, real_t t, real_t* U, real_t* U_dot, real_t* F)
{
  START_FUNCTION_TIMER();
  memcpy(integ->U_with_ghosts, U, sizeof(real_t) * integ->num_local_values);
  memcpy(integ->U_dot_with_ghosts, U_dot, sizeof(real_t) * integ->num_local_values);
  integ->F(integ->context, t, integ->U_with_ghosts, integ->U_dot_with_ghosts, F);
  STOP_FUNCTION_TIMER();
}

void dae_solver_set_max_dt(dae_solver_t* integ, real_t max_dt)
{
  ASSERT(max_dt > 0);
  integ->max_dt = max_dt;
  IDASetMaxStep(integ->ida, max_dt);
}

void dae_solver_set_stop_time(dae_solver_t* integ, real_t stop_time)
{
  integ->stop_time = stop_time;
  IDASetStopTime(integ->ida, stop_time);
}

bool dae_solver_step(dae_solver_t* integ, real_t max_dt, real_t* t, real_t* U, real_t* U_dot)
{
  START_FUNCTION_TIMER();
  int status = IDA_SUCCESS;

  // if we haven't been initialized, we need to copy in the data.
  if (!integ->initialized)
    dae_solver_reset(integ, *t, U, U_dot, DAE_IC_ASSUME_CONSISTENT);

  // If *t + max_dt is less than the time to which we've already integrated, 
  // we don't need to integrate; we only need to interpolate backward.
  real_t t2 = *t + max_dt;
  if (t2 > integ->t)
  {
    // Integrate to at least t -> t + max_dt.
    status = IDASolve(integ->ida, t2, &integ->t, integ->U, integ->U_dot, IDA_ONE_STEP);
    if ((status != IDA_SUCCESS) && (status != IDA_TSTOP_RETURN))
    {
      integ->status_message = get_status_message(status, integ->t);
      STOP_FUNCTION_TIMER();
      return false;
    }

    if ((t2 - *t) < (integ->t - *t))
      log_detail("dae_solver: took internal step dt = %g", integ->t - *t);
  }

  // If we integrated past t2, interpolate to t2.
  if (integ->t > t2)
  {
    status = IDAGetDky(integ->ida, t2, 0, integ->U);
    status = IDAGetDky(integ->ida, t2, 1, integ->U_dot);
    *t = t2;
  }
  else
    *t = integ->t;
  
  // Clear the present status.
  if (integ->status_message != NULL)
  {
    polymec_free(integ->status_message);
    integ->status_message = NULL;
  }

  // Did it work?
  if ((status == IDA_SUCCESS) || (status == IDA_TSTOP_RETURN))
  {
    // Copy out the solution.
    memcpy(U, NV_DATA(integ->U), sizeof(real_t) * integ->num_local_values); 
    memcpy(U_dot, NV_DATA(integ->U_dot), sizeof(real_t) * integ->num_local_values); 
    STOP_FUNCTION_TIMER();
    return true;
  }
  else
  {
    integ->status_message = get_status_message(status, integ->t);
    STOP_FUNCTION_TIMER();
    return false;
  }
}

void dae_solver_reset(dae_solver_t* integ, 
                      real_t t, 
                      real_t* U, 
                      real_t* U_dot,
                      dae_ic_correction_t ic_correction)
{
  // Reset the preconditioner.
  if (integ->precond != NULL)
    newton_pc_reset(integ->precond, t);

  // Reset the solver itself.
  integ->t = t;
  memcpy(NV_DATA(integ->U), U, sizeof(real_t) * integ->num_local_values); 
  memcpy(NV_DATA(integ->U_dot), U_dot, sizeof(real_t) * integ->num_local_values); 
  IDAReInit(integ->ida, integ->t, integ->U, integ->U_dot);
  integ->initialized = true;

  if (ic_correction != DAE_IC_ASSUME_CONSISTENT)
  {
    int err;
    if (ic_correction == DAE_IC_CORRECT_DERIVATIVES)
    {
      // Correct the initial conditions for U_dot given U.
      log_detail("dae_solver: computing initial condition derivatives...");
      err = IDACalcIC(integ->ida, IDA_YA_YDP_INIT, integ->max_dt);
    }
    else // if (ic_correction == DAE_IC_ASSUME_QUASISTATIC)
    {
      // Correct the initial conditions for U given U_dot = 0.
      log_detail("dae_solver: assuming quasistatic initial conditions...");
      err = IDACalcIC(integ->ida, IDA_Y_INIT, integ->max_dt);
    }
    if (err != IDA_SUCCESS)
    {
      log_detail("dae_solver: could not correct initial conditions:");
      if (err == IDA_LSETUP_FAIL)
        log_detail("  linear solver setup failed unrecoverably.");
      else if (err == IDA_LINIT_FAIL)
        log_detail("  linear solver initialization failed.");
      else if (err == IDA_LSOLVE_FAIL)
        log_detail("  linear solve failed.");
      else if (err == IDA_BAD_EWT)
        log_detail("  zero component of error weight function due to IC correction.");
      else if (err == IDA_FIRST_RES_FAIL)
        log_detail("  first call to residual failed, causing IC correction to fail.");
      else if (err == IDA_RES_FAIL)
        log_detail("  call to residual failed unrecoverably.");
      else if (err == IDA_CONSTR_FAIL)
        log_detail("  IC corrections were unable to meet constraints.");
      else if (err == IDA_LINESEARCH_FAIL)
        log_detail("  IC correction failed in line search.");
      else
        log_detail("  IC correction failed to converge.");
    }
    IDAGetConsistentIC(integ->ida, integ->U, integ->U_dot);
    memcpy(U, NV_DATA(integ->U), sizeof(real_t) * integ->num_local_values); 
    memcpy(U_dot, NV_DATA(integ->U_dot), sizeof(real_t) * integ->num_local_values); 
  }
  else
    log_detail("dae_solver: Assuming consistent initial conditions...");

  // Write out debugging info.
  if (log_level() == LOG_DEBUG)
  {
    long num_backtracks;
    IDAGetNumBacktrackOps(integ->ida, &num_backtracks);
    log_debug("dae_solver: backtracked %ld times in IC correction.", num_backtracks);
  }
}

void dae_solver_get_diagnostics(dae_solver_t* integ, 
                                dae_solver_diagnostics_t* diagnostics)
{
  diagnostics->status_message = integ->status_message; // borrowed!
  IDAGetNumSteps(integ->ida, &diagnostics->num_steps);
  IDAGetLastOrder(integ->ida, &diagnostics->order_of_last_step);
  IDAGetActualInitStep(integ->ida, &diagnostics->initial_step_size);
  IDAGetLastStep(integ->ida, &diagnostics->last_step_size);
  IDAGetNumResEvals(integ->ida, &diagnostics->num_residual_evaluations);
  IDAGetNumLinSolvSetups(integ->ida, &diagnostics->num_linear_solve_setups);
  IDAGetNumErrTestFails(integ->ida, &diagnostics->num_error_test_failures);
  IDAGetNumNonlinSolvIters(integ->ida, &diagnostics->num_nonlinear_solve_iterations);
  IDAGetNumNonlinSolvConvFails(integ->ida, &diagnostics->num_nonlinear_solve_convergence_failures);
  if (integ->solve_func == NULL) // JFNK mode
  {
    IDASpilsGetNumLinIters(integ->ida, &diagnostics->num_linear_solve_iterations);
    IDASpilsGetNumPrecEvals(integ->ida, &diagnostics->num_preconditioner_evaluations);
    IDASpilsGetNumPrecSolves(integ->ida, &diagnostics->num_preconditioner_solves);
    IDASpilsGetNumConvFails(integ->ida, &diagnostics->num_linear_solve_convergence_failures);
  }
  else
  {
    diagnostics->num_linear_solve_iterations = (long int)integ->num_linear_iterations;
    diagnostics->num_linear_solve_convergence_failures = (long int)integ->num_linear_conv_failures;
    diagnostics->num_preconditioner_evaluations = -1;
    diagnostics->num_preconditioner_solves = -1;
  }
}

void dae_solver_diagnostics_fprintf(dae_solver_diagnostics_t* diagnostics, 
                                    FILE* stream)
{
  if (stream == NULL) return;
  fprintf(stream, "dae solver diagnostics:\n");
  if (diagnostics->status_message != NULL)
    fprintf(stream, "  Status: %s\n", diagnostics->status_message);
  fprintf(stream, "  Num steps: %d\n", (int)diagnostics->num_steps);
  fprintf(stream, "  Order of last step: %d\n", diagnostics->order_of_last_step);
  fprintf(stream, "  Initial step size: %g\n", diagnostics->initial_step_size);
  fprintf(stream, "  Last step size: %g\n", diagnostics->last_step_size);
  fprintf(stream, "  Num residual evaluations: %d\n", (int)diagnostics->num_residual_evaluations);
  fprintf(stream, "  Num linear solve setups: %d\n", (int)diagnostics->num_linear_solve_setups);
  fprintf(stream, "  Num linear solve iterations: %d\n", (int)diagnostics->num_linear_solve_iterations);
  fprintf(stream, "  Num linear solve convergence failures: %d\n", (int)diagnostics->num_linear_solve_convergence_failures);
  fprintf(stream, "  Num error test failures: %d\n", (int)diagnostics->num_error_test_failures);
  fprintf(stream, "  Num nonlinear solve iterations: %d\n", (int)diagnostics->num_nonlinear_solve_iterations);
  fprintf(stream, "  Num nonlinear solve convergence failures: %d\n", (int)diagnostics->num_nonlinear_solve_convergence_failures);
  if (diagnostics->num_preconditioner_evaluations != -1) // JFNK mode
    fprintf(stream, "  Num preconditioner evaluations: %d\n", (int)diagnostics->num_preconditioner_evaluations);
  if (diagnostics->num_preconditioner_solves != -1) // JFNK mode
    fprintf(stream, "  Num preconditioner solves: %d\n", (int)diagnostics->num_preconditioner_solves);
}

//------------------------------------------------------------------------
//                  Inexact Newton-Krylov solver stuff
//------------------------------------------------------------------------

typedef struct
{
  MPI_Comm comm;

  // Behavior.
  void* context;
  int (*F_func)(void* context, real_t t, real_t* U, real_t* U_dot, real_t* F);
  int (*J_func)(void* context, real_t t, real_t* U, real_t alpha, real_t* U_dot, real_t* F, krylov_matrix_t* J);
  void (*dtor)(void* context);

  // Linear system stuff.
  krylov_factory_t* factory;
  krylov_solver_t* solver;
  krylov_pc_t* pc;
  krylov_matrix_t* J; // Newton matrix J = dF/dU + alpha * dF/d(U_dot)
  krylov_vector_t* X; // Solution vector.
  krylov_vector_t* B; // Right-hand side vector.
  krylov_vector_t* W; // Weights vector.

  // Metadata.
  matrix_sparsity_t* sparsity;
  int block_size;
} ink_dae_t;

static int ink_F(void* context, real_t t, real_t* U, real_t* U_dot, real_t* F)
{
  ink_dae_t* ink = context;
  return ink->F_func(ink->context, t, U, U_dot, F);
}

static int ink_reset(void* context, real_t t, real_t* U, real_t* U_dot)
{
  START_FUNCTION_TIMER();
  ink_dae_t* ink = context;
  
  // Free any resources currently in use.
  if (ink->J != NULL)
    krylov_matrix_free(ink->J);
  if (ink->block_size > 1)
    ink->J = krylov_factory_block_matrix(ink->factory, ink->sparsity, ink->block_size);
  else
    ink->J = krylov_factory_matrix(ink->factory, ink->sparsity);
  if (ink->X != NULL)
    krylov_vector_free(ink->X);
  if (ink->B != NULL)
    krylov_vector_free(ink->B);

  // Allocate resources.
  index_t* row_dist = matrix_sparsity_row_distribution(ink->sparsity);
  ink->X = krylov_factory_vector(ink->factory, ink->comm, row_dist);
  ink->B = krylov_factory_vector(ink->factory, ink->comm, row_dist);
  ink->W = krylov_factory_vector(ink->factory, ink->comm, row_dist);

  STOP_FUNCTION_TIMER();
  return 0;
}

static int ink_setup(void* context, 
                     real_t alpha, 
                     int step,
                     real_t t, 
                     real_t* U_pred, 
                     real_t* U_dot_pred, 
                     real_t* F_pred, 
                     real_t* work1, real_t* work2, real_t* work3)
{
  START_FUNCTION_TIMER();
  ink_dae_t* ink = context;

  log_debug("ink_dae_solver: Calculating J = dF/dU + %g * dF/d(U_dot).", alpha);
  int status = ink->J_func(ink->context, t, U_pred, alpha, U_dot_pred, F_pred, ink->J);
  if (status != 0)
    return status;

  // Use this matrix as the operator in our solver.
  krylov_solver_set_operator(ink->solver, ink->J);

  STOP_FUNCTION_TIMER();
  return 0;
}

static int ink_solve(void* context, 
                     real_t t, 
                     real_t* U,
                     real_t* U_dot,
                     real_t* F,
                     real_t* W, 
                     real_t res_norm_tol,
                     real_t* B,
                     int* num_iters) 
{
  START_FUNCTION_TIMER();
  ink_dae_t* ink = context;

  // Copy RHS data from B into ink->B.
  krylov_vector_copy_in(ink->B, B);

  // Copy weights into ink->W.
  krylov_vector_copy_in(ink->W, W);

  // If the WRMS norm of B is less than our tolerance, return X = 0.
  real_t B_norm = krylov_vector_norm(ink->B, 2);
  if (B_norm < res_norm_tol)
  {
    log_debug("ink_dae_solver: ||B||_2 < tolerance (%g < %g), so X -> 0.", B_norm, res_norm_tol); 
    krylov_vector_zero(ink->X);
    krylov_vector_copy_out(ink->X, B);
    STOP_FUNCTION_TIMER();
    return 0;
  }

  // Set the tolerances on the residual norm so that the relative 
  // tolerance is as restrictive as the absolute tolerance.
  real_t rel_tol = res_norm_tol / B_norm;
  real_t div_tol = 10.0;
  krylov_solver_set_tolerances(ink->solver, rel_tol, res_norm_tol, div_tol);

  // Solve A*X = B.
  real_t res_norm;
  bool solved = krylov_solver_solve_scaled(ink->solver, ink->B, ink->W, ink->W, 
                                           ink->X, &res_norm, num_iters);

  if (solved)
  {
    log_debug("ink_dae_solver: Solved A*X = B (||R|| == %g after %d iters).", res_norm, num_iters);

    // Copy solution data from ink->X into B.
    krylov_vector_copy_out(ink->X, B);
    STOP_FUNCTION_TIMER();
    return 0;
  }
  else
  {
    log_debug("ink_dae_solver: Solution to A*X = B did not converge.");
    STOP_FUNCTION_TIMER();
    return 1;
  }
}

static void ink_dtor(void* context)
{
  ink_dae_t* ink = context;
  matrix_sparsity_free(ink->sparsity);
  if (ink->X != NULL)
    krylov_vector_free(ink->X);
  if (ink->B != NULL)
    krylov_vector_free(ink->B);
  if (ink->W != NULL)
    krylov_vector_free(ink->W);
  if (ink->J != NULL)
    krylov_matrix_free(ink->J);
  if (ink->pc != NULL)
    krylov_pc_free(ink->pc);
  krylov_solver_free(ink->solver);
  krylov_factory_free(ink->factory);
  if ((ink->context != NULL) && (ink->dtor != NULL))
    ink->dtor(ink->context);
}

dae_solver_t* ink_dae_solver_new(int order, 
                                 MPI_Comm comm,
                                 dae_equation_t* equation_types,
                                 dae_constraint_t* constraints,
                                 krylov_factory_t* factory,
                                 matrix_sparsity_t* J_sparsity,
                                 void* context, 
                                 int (*F_func)(void* context, real_t t, real_t* U, real_t* U_dot, real_t* F),
                                 int (*J_func)(void* context, real_t t, real_t* U, real_t alpha, real_t* U_dot, real_t* F, krylov_matrix_t* J),
                                 void (*dtor)(void* context))
{
  ink_dae_t* ink = polymec_malloc(sizeof(ink_dae_t));
  ink->comm = comm;
  ink->context = context;
  ink->F_func = F_func;
  ink->J_func = J_func;
  ink->dtor = dtor;
  ink->factory = factory;
  ink->sparsity = J_sparsity;
  ink->solver = krylov_factory_gmres_solver(ink->factory, comm, 30);
  ink->pc = NULL;
  ink->J = NULL;
  ink->B = NULL;
  ink->X = NULL;
  ink->W = NULL;
  ink->block_size = 1;

  char name[1024];
  snprintf(name, 1024, "INK Differential Algebraic Equation solver (order %d)", order);
  int num_local_values = (int)(matrix_sparsity_num_local_rows(J_sparsity));
  dae_solver_t* I = dae_solver_new(name, order, comm, 
                                   equation_types, constraints,
                                   num_local_values, 0,
                                   ink, ink_F, ink_reset, 
                                   ink_setup, ink_solve, ink_dtor);

  // Set default tolerances.
  // relative error of 1e-4 means errors are controlled to 0.01%.
  // absolute error is set to 1 because it's completely problem dependent.
  dae_solver_set_tolerances(I, 1e-4, 1.0);

  return I;
}

void ink_dae_solver_use_pcg(dae_solver_t* ink_dae_integ)
{
  ink_dae_t* ink = dae_solver_context(ink_dae_integ);
  if (ink->solver != NULL)
    krylov_solver_free(ink->solver);
  ink->solver = krylov_factory_pcg_solver(ink->factory, ink->comm);
}

void ink_dae_solver_use_gmres(dae_solver_t* ink_dae_integ,
                              int max_krylov_dim)
{
  ink_dae_t* ink = dae_solver_context(ink_dae_integ);
  if (ink->solver != NULL)
    krylov_solver_free(ink->solver);
  ink->solver = krylov_factory_gmres_solver(ink->factory, ink->comm, max_krylov_dim);
}

void ink_dae_solver_use_bicgstab(dae_solver_t* ink_dae_integ)
{
  ink_dae_t* ink = dae_solver_context(ink_dae_integ);
  if (ink->solver != NULL)
    krylov_solver_free(ink->solver);
  ink->solver = krylov_factory_bicgstab_solver(ink->factory, ink->comm);
}

void ink_dae_solver_use_special(dae_solver_t* ink_dae_integ,
                                const char* solver_name,
                                string_string_unordered_map_t* options)
{
  ink_dae_t* ink = dae_solver_context(ink_dae_integ);
  if (ink->solver != NULL)
    krylov_solver_free(ink->solver);
  ink->solver = krylov_factory_special_solver(ink->factory, ink->comm,
                                              solver_name, options);
}

void ink_dae_solver_set_pc(dae_solver_t* ink_dae_integ,
                           const char* pc_name, 
                           string_string_unordered_map_t* options)
{
  ink_dae_t* ink = dae_solver_context(ink_dae_integ);
  if (ink->pc != NULL)
    krylov_pc_free(ink->pc);
  ink->pc = krylov_factory_preconditioner(ink->factory, ink->comm, pc_name, options);
}

void ink_dae_solver_set_block_size(dae_solver_t* ink_dae_integ,
                                   int block_size)
{
  ink_dae_t* ink = dae_solver_context(ink_dae_integ);
  ink->block_size = block_size;
}

void* ink_dae_solver_context(dae_solver_t* ink_dae_integ)
{
  ink_dae_t* ink = dae_solver_context(ink_dae_integ);
  return ink->context;
}

