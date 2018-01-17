// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <float.h>
#include "core/sundials_helpers.h"
#include "core/timer.h"
#include "solvers/newton_solver.h"

// We use KINSOL for doing the matrix-free nonlinear solve.
#include "kinsol/kinsol.h"
#include "kinsol/kinsol_impl.h"
#include "kinsol/kinsol_spgmr.h"
#include "kinsol/kinsol_spfgmr.h"
#include "kinsol/kinsol_spbcgs.h"
#include "kinsol/kinsol_sptfqmr.h"

struct newton_solver_t 
{
  // Parallel stuff.
  int rank, nprocs;
  MPI_Comm comm;
  int num_local_values, num_remote_values;

  // Behavior/state stuff.
  void* context;
  int (*F_func)(void* context, real_t t, real_t* U, real_t* F);
  void (*dtor)(void* context);

  // JFNK stuff.
  int (*Jv_func)(void* context, bool new_U, real_t t, real_t* U, real_t* v, real_t* Jv);
  jfnk_newton_t solver_type;
  newton_pc_t* precond;
  int max_krylov_dim, max_restarts;

  // Generized adaptor stuff.
  int (*reset_func)(void* context);
  int (*newton_setup_func)(void* context, 
                           newton_solver_strategy_t strategy,
                           real_t t, 
                           real_t* U,
                           real_t* F);
  int (*picard_setup_func)(void* context, 
                           real_t t, 
                           real_t* U,
                           real_t* F);
  int (*solve_func)(void* context, 
                    real_t* DF, 
                    real_t t, 
                    real_t* U,
                    real_t* F,
                    real_t* B,
                    real_t res_norm_tol,
                    real_t* p,
                    real_t* Jp_norm, 
                    real_t* F_o_Jp,
                    int* num_iters);
  int num_linear_iterations, num_linear_conv_failures;

  // KINSol data structures.
  void* kinsol;
  int strategy; // Global strategy.
  N_Vector U, U_scale, F_scale, constraints; // Stores solution vector and scaling vectors.
  real_t* U_with_ghosts;
  char* status_message; // status of most recent integration.

  // Current simulation time.
  real_t t;
};

// This function wraps around the user-supplied evaluation function.
static int evaluate_F(N_Vector U, N_Vector F, void* context)
{
  newton_solver_t* solver = context;
  real_t* UU = NV_DATA(U);
  real_t* FF = NV_DATA(F);

  // Evaluate the residual using a solution vector with ghosts.
  memcpy(solver->U_with_ghosts, UU, sizeof(real_t) * solver->num_local_values);
  return solver->F_func(solver->context, solver->t, solver->U_with_ghosts, FF);
}

// This function sets up the preconditioner data within the solver.
static int set_up_preconditioner(N_Vector U, N_Vector U_scale, 
                                 N_Vector F, N_Vector F_scale,
                                 void* context, 
                                 N_Vector work1, N_Vector work2)
{
  newton_solver_t* solver = context;
  real_t t = solver->t;
  newton_pc_setup(solver->precond, 0.0, 1.0, 0.0, t, NV_DATA(U), NULL);

  return 0;
}

// This function solves the preconditioner equation. On input, the vector r 
// contains the right-hand side of the preconditioner system, and on output 
// it contains the solution to the system.
static int solve_preconditioner_system(N_Vector U, N_Vector U_scale,
                                       N_Vector F, N_Vector F_scale,
                                       N_Vector r, void* context,
                                       N_Vector work)
{
  newton_solver_t* solver = context;

  // FIXME: Apply scaling if needed.

  if (newton_pc_solve(solver->precond, solver->t, NV_DATA(U), NULL,
                      NV_DATA(r), NV_DATA(work)))
  {
    // Copy the solution to r.
    memcpy(NV_DATA(r), NV_DATA(work), sizeof(real_t) * NV_LOCLENGTH(r));
    return 0;
  }
  else 
  {
    // Recoverable error.
    log_debug("newton_solver: preconditioner solve failed.");
    return 1; 
  }
}

static int newton_linit(KINMem kin_mem)
{
  newton_solver_t* newton = kin_mem->kin_user_data;
  newton->num_linear_iterations = 0;
  newton->num_linear_conv_failures = 0;
  return newton->reset_func(newton->context);
}

static int newton_lsetup(KINMem kin_mem)
{
  newton_solver_t* newton = kin_mem->kin_user_data;
  real_t t = newton->t;
  newton_solver_strategy_t strategy = NEWTON_FULL_STEP;
  if (newton->strategy == KIN_LINESEARCH) 
    strategy = NEWTON_LINE_SEARCH;
  int result = newton->newton_setup_func(newton->context, strategy, t, 
                                         NV_DATA(kin_mem->kin_uu), 
                                         NV_DATA(kin_mem->kin_fval));
  return result;
}

static int picard_lsetup(KINMem kin_mem)
{
  newton_solver_t* newton = kin_mem->kin_user_data;
  real_t t = newton->t;
  int result = newton->picard_setup_func(newton->context, t, 
                                         NV_DATA(kin_mem->kin_uu), 
                                         NV_DATA(kin_mem->kin_fval));
  return result;
}

static int newton_lsolve(KINMem kin_mem, 
                         N_Vector U, 
                         N_Vector B, 
                         real_t* sJpnorm,
                         real_t* sFdotJp)
{
  newton_solver_t* newton = kin_mem->kin_user_data;
  real_t t = newton->t;
  real_t* sJpnorm_ptr = NULL;
  real_t* sFdotJp_ptr = NULL;
  // We can shave a bit of work off by not computing norms unnecessarily.
  int strategy = newton->strategy;
  int eta_choice = kin_mem->kin_etaflag;
  if (eta_choice == KIN_ETACHOICE1)
  {
    sJpnorm_ptr = sJpnorm;
    sFdotJp_ptr = sFdotJp;
  }
  else if (strategy == KIN_LINESEARCH) 
    sFdotJp_ptr = sFdotJp;
  int num_iters = 0;
  int result =  newton->solve_func(newton->context, 
                                   NV_DATA(kin_mem->kin_fscale), t,
                                   NV_DATA(kin_mem->kin_uu), NV_DATA(kin_mem->kin_fval),
                                   NV_DATA(B), kin_mem->kin_eps, 
                                   NV_DATA(U), sJpnorm_ptr, sFdotJp_ptr, &num_iters);
  newton->num_linear_iterations += num_iters;
  if (result != 0)
    ++(newton->num_linear_conv_failures);
  return result;
}

static int newton_lfree(KINMem kin_mem)
{
  return 0;
}

static newton_solver_t* create_newton_solver(MPI_Comm comm,
                                             int num_local_values,
                                             int num_remote_values,
                                             void* context,
                                             int (*F_func)(void* context, real_t t, real_t* U, real_t* F),
                                             void (*dtor)(void* context))
{
  ASSERT(num_local_values > 0);
  ASSERT(num_remote_values >= 0);
  ASSERT(F_func != NULL);

  newton_solver_t* solver = polymec_malloc(sizeof(newton_solver_t));
  solver->context = context;
  solver->comm = comm;
  solver->F_func = F_func;
  solver->Jv_func = NULL;
  solver->reset_func = NULL;
  solver->newton_setup_func = NULL;
  solver->picard_setup_func = NULL;
  solver->solve_func = NULL;
  solver->dtor = dtor;
  solver->precond = NULL;
  solver->num_local_values = num_local_values;
  solver->num_remote_values = num_remote_values;
  solver->max_krylov_dim = -1;
  solver->max_restarts = -1;
  solver->strategy = KIN_NONE;

  // Set up KINSol and accessories.
  solver->kinsol = KINCreate();
  KINSetUserData(solver->kinsol, solver);
  solver->U = N_VNew(solver->comm, num_local_values);
  solver->U_with_ghosts = polymec_malloc(sizeof(real_t) * (solver->num_local_values + solver->num_remote_values));
  solver->U_scale = N_VNew(solver->comm, num_local_values);
  solver->F_scale = N_VNew(solver->comm, num_local_values);
  solver->constraints = N_VNew(solver->comm, num_local_values);
  solver->status_message = NULL;
  solver->t = 0.0;

  KINInit(solver->kinsol, evaluate_F, solver->U);

  // Set trivial scaling by default.
  newton_solver_set_U_scale(solver, NULL);
  newton_solver_set_F_scale(solver, NULL);

  // Enable debugging diagnostics if logging permits.
  FILE* info_stream = log_stream(LOG_DEBUG);
  if (info_stream != NULL)
  {
    KINSetPrintLevel(solver->kinsol, 3);
    KINSetInfoFile(solver->kinsol, info_stream);
  }
  else
  {
    KINSetPrintLevel(solver->kinsol, 0);
    KINSetInfoFile(solver->kinsol, NULL);
  }

  return solver;
}

newton_solver_t* newton_solver_new(MPI_Comm comm,
                                   int num_local_values,
                                   int num_remote_values,
                                   void* context,
                                   int (*F_func)(void* context, real_t t, real_t* U, real_t* F),
                                   int (*reset_func)(void* context),
                                   int (*setup_func)(void* context, 
                                                     newton_solver_strategy_t strategy,
                                                     real_t t,
                                                     real_t* U,
                                                     real_t* F),
                                   int (*solve_func)(void* context, 
                                                     real_t* DF, 
                                                     real_t t, 
                                                     real_t* U,
                                                     real_t* F,
                                                     real_t* B,
                                                     real_t res_norm_tol,
                                                     real_t* p,
                                                     real_t* Jp_norm, 
                                                     real_t* F_o_Jp,
                                                     int* num_iters), 
                                   void (*dtor)(void* context),
                                   newton_solver_strategy_t strategy)
{
  ASSERT(reset_func != NULL);
  ASSERT(setup_func != NULL);
  ASSERT(solve_func != NULL);

  newton_solver_t* solver = create_newton_solver(comm,
                                                 num_local_values,
                                                 num_remote_values,
                                                 context,
                                                 F_func,
                                                 dtor);

  // Fill in the holes.
  solver->reset_func = reset_func;
  solver->newton_setup_func = setup_func;
  solver->solve_func = solve_func;
  solver->strategy = (strategy == NEWTON_FULL_STEP) ? KIN_NONE : KIN_LINESEARCH;

  // Set up our generalized solver.
  KINMem kin_mem = solver->kinsol;
  kin_mem->kin_linit = newton_linit;
  kin_mem->kin_lsetup = newton_lsetup;
  kin_mem->kin_lsolve = newton_lsolve;
  kin_mem->kin_lfree = newton_lfree;
  kin_mem->kin_setupNonNull = 1;  // need this for lsetup to be called
  kin_mem->kin_inexact_ls = 1;    // need this for iterative solvers

  return solver;
}

newton_solver_t* picard_newton_solver_new(MPI_Comm comm,
                                          int num_local_values,
                                          int num_remote_values,
                                          void* context,
                                          int (*F_func)(void* context, real_t t, real_t* U, real_t* F),
                                          int (*reset_func)(void* context),
                                          int (*setup_func)(void* context, 
                                                            real_t t,
                                                            real_t* U,
                                                            real_t* F),
                                          int (*solve_func)(void* context, 
                                                            real_t* DF, 
                                                            real_t t, 
                                                            real_t* U,
                                                            real_t* F,
                                                            real_t* B,
                                                            real_t res_norm_tol,
                                                            real_t* p,
                                                            real_t* Lp_norm, 
                                                            real_t* F_o_Lp,
                                                            int* num_iters), 
                                          void (*dtor)(void* context),
                                          int num_residuals)
{
  ASSERT(reset_func != NULL);
  ASSERT(setup_func != NULL);
  ASSERT(solve_func != NULL);

  newton_solver_t* solver = create_newton_solver(comm,
                                                 num_local_values,
                                                 num_remote_values,
                                                 context,
                                                 F_func,
                                                 dtor);

  // Fill in the holes.
  solver->reset_func = reset_func;
  solver->picard_setup_func = setup_func;
  solver->solve_func = solve_func;
  solver->strategy = KIN_PICARD;
  KINSetMAA(solver->kinsol, num_residuals);

  // Set up our generalized solver.
  KINMem kin_mem = solver->kinsol;
  kin_mem->kin_linit = newton_linit;
  kin_mem->kin_lsetup = picard_lsetup;
  kin_mem->kin_lsolve = newton_lsolve;
  kin_mem->kin_lfree = newton_lfree;
  kin_mem->kin_setupNonNull = 1;  // need this for lsetup to be called
  kin_mem->kin_inexact_ls = 1;    // need this for iterative solvers

  return solver;
}

newton_solver_t* fixed_point_newton_solver_new(MPI_Comm comm,
                                               int num_local_values,
                                               int num_remote_values,
                                               void* context,
                                               int (*G_func)(void* context, real_t t, real_t* U, real_t* G),
                                               void (*dtor)(void* context),
                                               int num_residuals)
{
  newton_solver_t* solver = create_newton_solver(comm,
                                                 num_local_values,
                                                 num_remote_values,
                                                 context,
                                                 G_func,
                                                 dtor);
  solver->strategy = KIN_FP;
  KINSetMAA(solver->kinsol, num_residuals);

  // Do we need this?
  KINSpbcg(solver->kinsol, 30);

  return solver;
}
       
// This is a KINSOL Jacobian-vector product function that wraps our own.
static int jfnk_Jv_func_wrapper(N_Vector v, N_Vector Jv, N_Vector U,
                                booleantype* new_U, void* context)
{
  newton_solver_t* solver = context;
  ASSERT(solver->Jv_func != NULL);
  int result = solver->Jv_func(solver->context, (*new_U != 0), solver->t, 
                               NV_DATA(U), NV_DATA(v), NV_DATA(Jv));
  if ((result == 0) && (*new_U == 0))
    *new_U = 1;
  return result;
}

newton_solver_t* jfnk_newton_solver_new(MPI_Comm comm,
                                        int num_local_values,
                                        int num_remote_values,
                                        void* context,
                                        int (*F_func)(void* context, real_t t, real_t* U, real_t* F),
                                        int (*Jv_func)(void* context, bool new_U, real_t t, real_t* U, real_t* v, real_t* Jv),
                                        void (*dtor)(void* context),
                                        newton_solver_strategy_t strategy,
                                        newton_pc_t* precond,
                                        jfnk_newton_t solver_type,
                                        int max_krylov_dim, 
                                        int max_restarts)
{
  ASSERT(max_krylov_dim >= 3);
  ASSERT(((solver_type != NEWTON_GMRES) && (solver_type != NEWTON_FGMRES)) || 
         (max_restarts >= 0));
  ASSERT(precond != NULL);
  ASSERT(newton_pc_side(precond) == NEWTON_PC_LEFT);

  newton_solver_t* solver = create_newton_solver(comm,
                                                 num_local_values,
                                                 num_remote_values,
                                                 context,
                                                 F_func,
                                                 dtor);
  solver->Jv_func = Jv_func;
  solver->precond = precond;
  solver->solver_type = solver_type;
  solver->max_krylov_dim = max_krylov_dim;
  solver->max_restarts = max_restarts;
  solver->strategy = (strategy == NEWTON_FULL_STEP) ? KIN_NONE : KIN_LINESEARCH;

  // Select the particular type of Krylov method for the underlying linear solves.
  if (solver->solver_type == NEWTON_GMRES)
  {
    KINSpgmr(solver->kinsol, solver->max_krylov_dim); 
    KINSpilsSetMaxRestarts(solver->kinsol, solver->max_restarts);
  }
  else if (solver->solver_type == NEWTON_FGMRES)
  {
    KINSpfgmr(solver->kinsol, solver->max_krylov_dim); 
    KINSpilsSetMaxRestarts(solver->kinsol, solver->max_restarts);
  }
  else if (solver->solver_type == NEWTON_BICGSTAB)
    KINSpbcg(solver->kinsol, solver->max_krylov_dim);
  else
    KINSptfqmr(solver->kinsol, solver->max_krylov_dim);

  // Set up the Jacobian-vector product.
  if (Jv_func != NULL)
    KINSpilsSetJacTimesVecFn(solver->kinsol, jfnk_Jv_func_wrapper);
  else
    KINSpilsSetJacTimesVecFn(solver->kinsol, NULL);

  // Set up the preconditioner.
  if (solver->precond != NULL)
  {
    KINSpilsSetPreconditioner(solver->kinsol, set_up_preconditioner,
                              solve_preconditioner_system);
  }

  return solver;
}

newton_solver_t* picard_jfnk_newton_solver_new(MPI_Comm comm,
                                               int num_local_values,
                                               int num_remote_values,
                                               void* context,
                                               int (*F_func)(void* context, real_t t, real_t* U, real_t* F),
                                               int (*Lv_func)(void* context, bool new_U, real_t t, real_t* U, real_t* v, real_t* Jv),
                                               void (*dtor)(void* context),
                                               newton_pc_t* precond,
                                               jfnk_newton_t solver_type,
                                               int max_krylov_dim,
                                               int max_restarts,
                                               int num_residuals)
{
  newton_solver_t* solver = jfnk_newton_solver_new(comm, num_local_values, 
                                                   num_remote_values, context, 
                                                   F_func, Lv_func, dtor, NEWTON_FULL_STEP, 
                                                   precond, solver_type, 
                                                   max_krylov_dim, max_restarts);
  solver->strategy = KIN_PICARD;
  KINSetMAA(solver->kinsol, num_residuals);
  return solver;
}

void newton_solver_free(newton_solver_t* solver)
{
  // Kill the preconditioner stuff.
  if (solver->precond != NULL)
    newton_pc_free(solver->precond);

  // Kill the KINSol stuff.
  N_VDestroy(solver->U);
  N_VDestroy(solver->U_scale);
  N_VDestroy(solver->F_scale);
  N_VDestroy(solver->constraints);
  KINFree(&solver->kinsol);

  // Kill the rest.
  if ((solver->dtor != NULL) && (solver->context != NULL))
    solver->dtor(solver->context);

  // Kill the rest.
  if (solver->status_message != NULL)
    polymec_free(solver->status_message);
  polymec_free(solver->U_with_ghosts);
  polymec_free(solver);
}

void* newton_solver_context(newton_solver_t* solver)
{
  return solver->context;
}

int newton_solver_num_equations(newton_solver_t* solver)
{
  return solver->num_local_values;
}

void newton_solver_set_U_scale(newton_solver_t* solver, 
                               real_t* DU)
{
  if (DU != NULL)
    memcpy(NV_DATA(solver->U_scale), DU, sizeof(real_t) * solver->num_local_values);
  else
  {
    real_t* U_scale = NV_DATA(solver->U_scale);
    for (int i = 0; i < solver->num_local_values; ++i)
      U_scale[i] = 1.0;
  }
}

void newton_solver_set_F_scale(newton_solver_t* solver, 
                               real_t* DF)
{
  if (DF != NULL)
    memcpy(NV_DATA(solver->F_scale), DF, sizeof(real_t) * solver->num_local_values);
  else
  {
    real_t* F_scale = NV_DATA(solver->F_scale);
    for (int i = 0; i < solver->num_local_values; ++i)
      F_scale[i] = 1.0;
  }
}

void newton_solver_set_constraints(newton_solver_t* solver,
                                   real_t* constraints)
{
  if (constraints != NULL)
    memcpy(NV_DATA(solver->constraints), constraints, sizeof(real_t) * solver->num_local_values);
  else
  {
    real_t* C = NV_DATA(solver->constraints);
    for (int i = 0; i < solver->num_local_values; ++i)
      C[i] = 0.0;
  }
  KINSetConstraints(solver->kinsol, solver->constraints);
}

void newton_solver_set_tolerances(newton_solver_t* solver, real_t norm_tolerance, real_t step_tolerance)
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

void newton_solver_set_linear_solver_stopping_criteria(newton_solver_t* solver,
                                                       newton_solver_stopping_criteria_t criteria)
{
  if (criteria == NEWTON_EISENSTAT_WALKER1)
    KINSetEtaForm(solver->kinsol, KIN_ETACHOICE1);
  else if (criteria == NEWTON_EISENSTAT_WALKER2)
    KINSetEtaForm(solver->kinsol, KIN_ETACHOICE2);
  else // (criteria == NEWTON_CONSTANT_ETA)
    KINSetEtaForm(solver->kinsol, KIN_ETACONSTANT);
}

void newton_solver_set_constant_eta(newton_solver_t* solver, real_t eta)
{
  ASSERT(eta > 0.0);
  KINSetEtaConstValue(solver->kinsol, eta);
}

newton_pc_t* newton_solver_preconditioner(newton_solver_t* solver)
{
  return solver->precond;
}

void newton_solver_eval_residual(newton_solver_t* solver, 
                                 real_t t, 
                                 real_t* U, 
                                 real_t* R)
{
  solver->F_func(solver->context, t, U, R);
}

bool newton_solver_solve(newton_solver_t* solver,
                         real_t t,
                         real_t* U,
                         int* num_iterations)
{
  START_FUNCTION_TIMER();
  ASSERT(U != NULL);

  // Set the current time in the state.
  solver->t = t;

  // Copy the values in U to the internal solution vector.
  memcpy(NV_DATA(solver->U), U, sizeof(real_t) * solver->num_local_values);

  // Suspend the currently active floating point exceptions for now.
//  polymec_suspend_fpe_exceptions();

  // Solve.
  log_debug("newton_solver: solving...");
  int status = KINSol(solver->kinsol, solver->U, solver->strategy, 
                      solver->U_scale, solver->F_scale);

  // Clear the present status.
  if (solver->status_message != NULL)
  {
    polymec_free(solver->status_message);
    solver->status_message = NULL;
  }

  // Reinstate the floating point exceptions.
//  polymec_restore_fpe_exceptions();

  if ((status == KIN_SUCCESS) || (status == KIN_INITIAL_GUESS_OK))
  {
    // Get the number of iterations it took.
    long num_iters;
    KINGetNumNonlinSolvIters(solver->kinsol, &num_iters);
    *num_iterations = (int)num_iters;
    log_debug("newton_solver: solved after %d iterations.", *num_iterations);

    // Copy the data back into U.
    memcpy(U, NV_DATA(solver->U), sizeof(real_t) * solver->num_local_values);
    STOP_FUNCTION_TIMER();
    return true;
  }
  else
  {
    if (status == KIN_STEP_LT_STPTOL)
      solver->status_message = string_dup("Nonlinear solve stalled because scaled Newton step is too small.");
    else if (status == KIN_LINESEARCH_NONCONV)
      solver->status_message = string_dup("Line search could not sufficiently decrease the error of the iterate.");
    else if (status == KIN_MAXITER_REACHED)
      solver->status_message = string_dup("Maximum number of nonlinear iterations was reached.");
    else if (status == KIN_MXNEWT_5X_EXCEEDED)
      solver->status_message = string_dup("Maximum Newton step size was exceeded 5 times.");
    else if (status == KIN_LINESEARCH_BCFAIL)
      solver->status_message = string_dup("Line search could not satisfy beta condition.");
    else if (status == KIN_LINSOLV_NO_RECOVERY)
      solver->status_message = string_dup("Preconditioner solve encountered a recoverable error after update.");
    else if (status == KIN_LINIT_FAIL)
      solver->status_message = string_dup("Linear solve setup failed.");
    else if (status == KIN_LSETUP_FAIL)
      solver->status_message = string_dup("Preconditioner setup failed unrecoverably.");
    else if (status == KIN_LSOLVE_FAIL)
      solver->status_message = string_dup("Linear solve failed (or preconditioner solve failed unrecoverably).");
    else if (status == KIN_SYSFUNC_FAIL)
      solver->status_message = string_dup("Nonlinear function evaluation failed unrecoverably.");
    else if (status == KIN_FIRST_SYSFUNC_ERR)
      solver->status_message = string_dup("First nonlinear function evaluation failed recoverably.");
    else if (status == KIN_REPTD_SYSFUNC_ERR)
      solver->status_message = string_dup("Nonlinear function evaluation repeatedly failed (no recovery possible).");
  }

  // Failed!
  STOP_FUNCTION_TIMER();
  return false;
}
                                  
void newton_solver_reset(newton_solver_t* solver, real_t t)
{
  // Reset the preconditioner if it exists.
  if (solver->precond != NULL)
    newton_pc_reset(solver->precond, t);
}

void newton_solver_get_diagnostics(newton_solver_t* solver, 
                                   newton_solver_diagnostics_t* diagnostics)
{
  diagnostics->status_message = solver->status_message; // borrowed!
  KINGetNumFuncEvals(solver->kinsol, &diagnostics->num_function_evaluations);
  KINGetNumBetaCondFails(solver->kinsol, &diagnostics->num_beta_condition_failures);
  KINGetNumNonlinSolvIters(solver->kinsol, &diagnostics->num_nonlinear_iterations);
  KINGetNumBacktrackOps(solver->kinsol, &diagnostics->num_backtrack_operations);
  KINGetFuncNorm(solver->kinsol, &diagnostics->scaled_function_norm);
  KINGetStepLength(solver->kinsol, &diagnostics->scaled_newton_step_length);
  if (solver->solve_func == NULL) // JFNK mode
  {
    KINSpilsGetNumLinIters(solver->kinsol, &diagnostics->num_linear_solve_iterations);
    KINSpilsGetNumConvFails(solver->kinsol, &diagnostics->num_linear_solve_convergence_failures);
    KINSpilsGetNumPrecEvals(solver->kinsol, &diagnostics->num_preconditioner_evaluations);
    KINSpilsGetNumPrecSolves(solver->kinsol, &diagnostics->num_preconditioner_solves);
    KINSpilsGetNumJtimesEvals(solver->kinsol, &diagnostics->num_jacobian_vector_product_evaluations);
    KINSpilsGetNumFuncEvals(solver->kinsol, &diagnostics->num_difference_quotient_function_evaluations);
  }
  else
  {
    diagnostics->num_linear_solve_iterations = (long int)solver->num_linear_iterations;
    diagnostics->num_linear_solve_convergence_failures = (long int)solver->num_linear_conv_failures;
    diagnostics->num_preconditioner_evaluations = -1;
    diagnostics->num_preconditioner_solves = -1;
    diagnostics->num_jacobian_vector_product_evaluations = -1;
    diagnostics->num_difference_quotient_function_evaluations = -1;
  }
}

void newton_solver_diagnostics_fprintf(newton_solver_diagnostics_t* diagnostics, 
                                       FILE* stream)
{
  if (stream == NULL) return;
  fprintf(stream, "Nonlinear solver diagnostics:\n");
  if (diagnostics->status_message != NULL)
    fprintf(stream, "  Status: %s\n", diagnostics->status_message);
  fprintf(stream, "  Num function evaluations: %d\n", (int)diagnostics->num_function_evaluations);
  fprintf(stream, "  Num beta condition failures: %d\n", (int)diagnostics->num_beta_condition_failures);
  fprintf(stream, "  Num backtrack operations: %d\n", (int)diagnostics->num_backtrack_operations);
  fprintf(stream, "  Num nonlinear iterations: %d\n", (int)diagnostics->num_nonlinear_iterations);
  fprintf(stream, "  Scaled function norm: %g\n", (double)diagnostics->scaled_function_norm);
  fprintf(stream, "  Scaled Newton step length: %g\n", (double)diagnostics->scaled_newton_step_length);
  fprintf(stream, "  Num linear solve iterations: %d\n", (int)diagnostics->num_linear_solve_iterations);
  fprintf(stream, "  Num linear solve convergence failures: %d\n", (int)diagnostics->num_linear_solve_convergence_failures);
  if (diagnostics->num_preconditioner_evaluations != -1)
    fprintf(stream, "  Num preconditioner evaluations: %d\n", (int)diagnostics->num_preconditioner_evaluations);
  if (diagnostics->num_preconditioner_solves != -1)
    fprintf(stream, "  Num preconditioner solves: %d\n", (int)diagnostics->num_preconditioner_solves);
  if (diagnostics->num_jacobian_vector_product_evaluations != -1)
    fprintf(stream, "  Num Jacobian-vector product evaluations: %d\n", (int)diagnostics->num_jacobian_vector_product_evaluations);
  if (diagnostics->num_difference_quotient_function_evaluations != -1)
    fprintf(stream, "  Num difference quotient function evaluations: %d\n", (int)diagnostics->num_difference_quotient_function_evaluations);
}

//------------------------------------------------------------------------
//                      Inexact Newton-Krylov solver
//------------------------------------------------------------------------

typedef struct
{
  MPI_Comm comm;

  // Behavior.
  void* context;
  int (*F_func)(void* context, real_t t, real_t* U, real_t* F);
  int (*J_func)(void* context, real_t t, real_t* U, real_t* F, krylov_matrix_t* J);
  void (*dtor)(void* context);

  // Linear system stuff.
  krylov_factory_t* factory;
  krylov_solver_t* solver;
  krylov_pc_t* pc;
  krylov_matrix_t* J; // Jacobian matrix J.
  krylov_vector_t* X; // Increment vector.
  krylov_vector_t* B; // Right-hand side vector.
  krylov_vector_t* DF; // Scaling vector for F.
  krylov_vector_t* F; // Residual/system vector.

  // Metadata.
  matrix_sparsity_t* sparsity;
  int block_size;
} ink_newton_t;

static int ink_F(void* context, real_t t, real_t* U, real_t* F)
{
  ink_newton_t* ink = context;
  return ink->F_func(ink->context, t, U, F);
}

static int ink_reset(void* context)
{
  START_FUNCTION_TIMER();
  ink_newton_t* ink = context;
  
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
  if (ink->DF != NULL)
    krylov_vector_free(ink->DF);
  if (ink->F != NULL)
    krylov_vector_free(ink->F);

  // Allocate resources.
  index_t* row_dist = matrix_sparsity_row_distribution(ink->sparsity);
  ink->X = krylov_factory_vector(ink->factory, ink->comm, row_dist);
  ink->B = krylov_factory_vector(ink->factory, ink->comm, row_dist);
  ink->DF = krylov_factory_vector(ink->factory, ink->comm, row_dist);
  ink->F = krylov_factory_vector(ink->factory, ink->comm, row_dist);

  STOP_FUNCTION_TIMER();
  return 0;
}

static int ink_newton_setup(void* context, 
                            newton_solver_strategy_t strategy,
                            real_t t,
                            real_t* U,
                            real_t* F)
{
  START_FUNCTION_TIMER();
  ink_newton_t* ink = context;
  log_debug("ink_newton_solver: Calculating J = dF/dU.");

  // Compute the matrix.
  int status = ink->J_func(ink->context, t, U, F, ink->J);
  if (status != 0)
    return status;

  // Use this matrix as the operator in our solver.
  krylov_solver_set_operator(ink->solver, ink->J);

  STOP_FUNCTION_TIMER();
  return 0;
}

static int ink_picard_setup(void* context, 
                            real_t t,
                            real_t* U,
                            real_t* F)
{
  START_FUNCTION_TIMER();
  ink_newton_t* ink = context;
  log_debug("ink_picard_newton_solver: Calculating L.");

  // Compute the matrix.
  int status = ink->J_func(ink->context, t, U, F, ink->J);
  if (status != 0)
    return status;

  // Use this matrix as the operator in our solver.
  krylov_solver_set_operator(ink->solver, ink->J);

  STOP_FUNCTION_TIMER();
  return 0;
}

static int ink_solve(void* context, 
                     real_t* DF, 
                     real_t t, 
                     real_t* U,
                     real_t* F,
                     real_t* B,
                     real_t res_norm_tol,
                     real_t* p,
                     real_t* Jp_norm, 
                     real_t* F_o_Jp,
                     int* num_iters)
{
  START_FUNCTION_TIMER();
  ink_newton_t* ink = context;

  // Copy RHS data from B into ink->B.
  krylov_vector_copy_in(ink->B, B);

  // Copy the F scaling into ink->DF.
  krylov_vector_copy_in(ink->DF, DF);

  // Set the initial guess to 0.
  krylov_vector_zero(ink->X);

  // Set the tolerance on the residual norm.
  real_t rel_tol = 1e-8;
  real_t div_tol = 10.0;
  krylov_solver_set_tolerances(ink->solver, rel_tol, res_norm_tol, div_tol);

  // Solve (DF*A)*X = DF*B.
  real_t res_norm;
  bool solved = krylov_solver_solve_scaled(ink->solver, ink->B, ink->DF, ink->DF, 
                                           ink->X, &res_norm, num_iters);

  if (solved)
  {
    log_debug("ink_newton_solver: Solved A*X = B (||DF*(B-A*X)||_2 == %g", res_norm);
    log_debug("ink_newton_solver:                 after %d iterations).", num_iters);

    // Compute the norms, using ink->B as a workspace.
    if ((Jp_norm != NULL) || (F_o_Jp != NULL))
      krylov_matrix_matvec(ink->J, ink->X, false, ink->B); // J*p -> B.
    if (Jp_norm != NULL)
    {
      *Jp_norm = krylov_vector_w2_norm(ink->B, ink->DF);
      log_debug("ink_newton_solver: ||DF*J*p||_2 = %g", *Jp_norm);
    }
    if (F_o_Jp != NULL)
    {
      krylov_vector_diag_scale(ink->B, ink->DF);
      krylov_vector_diag_scale(ink->B, ink->DF);
      krylov_vector_copy_in(ink->F, F);
      *F_o_Jp = krylov_vector_dot(ink->F, ink->B);
      log_debug("ink_newton_solver: (DF*F) o (DF*J*p) = %g", *F_o_Jp);
    }

    // Copy solution data from ink->X into p.
    krylov_vector_copy_out(ink->X, p);
    STOP_FUNCTION_TIMER();
    return 0;
  }
  else
  {
    log_debug("ink_bdf_ode_solver: Solution to A*X = B did not converge.");
    STOP_FUNCTION_TIMER();
    return 1;
  }
}

static void ink_dtor(void* context)
{
  ink_newton_t* ink = context;
  matrix_sparsity_free(ink->sparsity);
  if (ink->F != NULL)
    krylov_vector_free(ink->F);
  if (ink->DF != NULL)
    krylov_vector_free(ink->DF);
  if (ink->X != NULL)
    krylov_vector_free(ink->X);
  if (ink->B != NULL)
    krylov_vector_free(ink->B);
  if (ink->J != NULL)
    krylov_matrix_free(ink->J);
  if (ink->pc != NULL)
    krylov_pc_free(ink->pc);
  krylov_solver_free(ink->solver);
  krylov_factory_free(ink->factory);
  if ((ink->context != NULL) && (ink->dtor != NULL))
    ink->dtor(ink->context);
}

newton_solver_t* ink_newton_solver_new(MPI_Comm comm,
                                       krylov_factory_t* factory,
                                       matrix_sparsity_t* J_sparsity,
                                       void* context, 
                                       int (*F_func)(void* context, real_t t, real_t* U, real_t* F),
                                       int (*J_func)(void* context, real_t t, real_t* U, real_t* F, krylov_matrix_t* J),
                                       void (*dtor)(void* context),
                                       newton_solver_strategy_t strategy)
{
  ink_newton_t* ink = polymec_malloc(sizeof(ink_newton_t));
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
  ink->DF = NULL;
  ink->F = NULL;
  ink->block_size = 1;

  int num_local_values = (int)(matrix_sparsity_num_local_rows(J_sparsity));
  newton_solver_t* N = newton_solver_new(comm, num_local_values, 0,
                                         ink, ink_F, ink_reset, 
                                         ink_newton_setup, ink_solve, ink_dtor,
                                         strategy);

  return N;
}

newton_solver_t* ink_picard_newton_solver_new(MPI_Comm comm,
                                              krylov_factory_t* factory,
                                              matrix_sparsity_t* L_sparsity,
                                              void* context, 
                                              int (*F_func)(void* context, real_t t, real_t* U, real_t* F),
                                              int (*L_func)(void* context, real_t t, real_t* U, real_t* F, krylov_matrix_t* J),
                                              void (*dtor)(void* context),
                                              int num_residuals)
{
  ink_newton_t* ink = polymec_malloc(sizeof(ink_newton_t));
  ink->comm = comm;
  ink->context = context;
  ink->F_func = F_func;
  ink->J_func = L_func;
  ink->dtor = dtor;
  ink->factory = factory;
  ink->sparsity = L_sparsity;
  ink->solver = krylov_factory_gmres_solver(ink->factory, comm, 30);
  ink->pc = NULL;
  ink->J = NULL;
  ink->B = NULL;
  ink->X = NULL;
  ink->DF = NULL;
  ink->F = NULL;
  ink->block_size = 1;

  int num_local_values = (int)(matrix_sparsity_num_local_rows(L_sparsity));
  newton_solver_t* N = picard_newton_solver_new(comm, num_local_values, 0,
                                                ink, ink_F, ink_reset, 
                                                ink_picard_setup, ink_solve, ink_dtor,
                                                num_residuals);

  return N;
}

