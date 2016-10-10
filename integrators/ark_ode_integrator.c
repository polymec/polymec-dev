// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/sundials_helpers.h"
#include "core/array.h"
#include "core/timer.h"
#include "integrators/ark_ode_integrator.h"
#include "integrators/newton_pc.h"
#include "arkode/arkode.h"

// JFNK stuff.
#include "arkode/arkode_spils.h"
#include "arkode/arkode_spgmr.h"
#include "arkode/arkode_spfgmr.h"
#include "arkode/arkode_spbcgs.h"
#include "arkode/arkode_pcg.h"
#include "arkode/arkode_sptfqmr.h"

// Stuff for generalized ARK integrators.
#include "arkode/arkode_impl.h"

struct ark_ode_observer_t 
{
  void* context;
  void (*fe_computed)(void* context, real_t t, real_t* U, real_t* fe);
  void (*fi_computed)(void* context, real_t t, real_t* U, real_t* fi);
  void (*Jy_computed)(void* context, real_t t, real_t* U, real_t* U_dot, real_t* y, real_t* Jy);
  void (*dtor)(void* context);
};

typedef struct
{
  MPI_Comm comm;
  int num_local_values, num_remote_values;
  void* context; 
  int (*fe)(void* context, real_t t, real_t* U, real_t* fe);
  int (*fi)(void* context, real_t t, real_t* U, real_t* fi);
  real_t (*stable_dt)(void* context, real_t t, real_t* U);
  void (*dtor)(void* context);

  // ARKode data structures.
  void* arkode;
  N_Vector U; 
  real_t* U_with_ghosts;
  real_t t;
  char* status_message; // status of most recent integration.

  bool first_step;

  // JFNK stuff.
  int max_krylov_dim;
  newton_pc_t* precond;
  int (*Jy)(void* context, real_t t, real_t* U, real_t* U_dot, real_t* y, real_t* temp, real_t* Jy);

  // Generalized adaptor stuff.
  real_t sqrtN;
  int (*reset_func)(void* context, real_t t, real_t* U);
  int (*setup_func)(void* context, 
                    ark_conv_status_t conv_status, 
                    real_t gamma, 
                    int step,
                    real_t t, 
                    real_t* U_pred, 
                    real_t* U_dot_pred, 
                    bool* J_current, 
                    real_t* work1, real_t* work2, real_t* work3);
  int (*solve_func)(void* context, 
                    real_t t, 
                    real_t* U,
                    real_t* U_dot,
                    real_t* W, 
                    real_t res_norm_tol,
                    real_t* B);

  // Error weight function.
  void (*compute_weights)(void* context, real_t* y, real_t* weights);
  real_t* error_weights;

  // Observers.
  ptr_array_t* observers;
} ark_ode_t;

static char* get_status_message(int status, real_t current_time)
{
  char* status_message = NULL;
  if (status == ARK_TOO_MUCH_WORK)
  {
    char err[1024];
    snprintf(err, 1024, "Integrator stopped at t = %g after maximum number of steps.", current_time);
    status_message = string_dup(err);
  }
  else if (status == ARK_TOO_MUCH_ACC)
    status_message = string_dup("Integrator could not achieve desired level of accuracy.");
  else if (status == ARK_ERR_FAILURE)
    status_message = string_dup("Integrator encountered too many error test failures.");
  else if (status == ARK_CONV_FAILURE)
    status_message = string_dup("Integrator encountered too many convergence test failures.");
  else if (status == ARK_LINIT_FAIL)
    status_message = string_dup("Integrator's linear solver failed to initialize.");
  else if (status == ARK_LSETUP_FAIL)
    status_message = string_dup("Integrator's linear solver setup failed.");
  else if (status == ARK_LSOLVE_FAIL)
    status_message = string_dup("Integrator's linear solver failed.");
  else if (status == ARK_RHSFUNC_FAIL)
    status_message = string_dup("Integrator's RHS function failed unrecoverably.");
  else if (status == ARK_FIRST_RHSFUNC_ERR)
    status_message = string_dup("Integrator's first call to RHS function failed.");
  else if (status == ARK_REPTD_RHSFUNC_ERR)
    status_message = string_dup("Integrator encountered too many recoverable RHS failures.");
  else if (status == ARK_UNREC_RHSFUNC_ERR)
    status_message = string_dup("Integrator failed to recover from a recoverable RHS failure.");
  else if (status == ARK_RTFUNC_FAIL)
    status_message = string_dup("Integrator encountered a failure in the rootfinding function.");
  return status_message;
}

// This function wraps around the user-supplied explicit RHS.
static int evaluate_fe(real_t t, N_Vector U, N_Vector U_dot, void* context)
{
  START_FUNCTION_TIMER();
  ark_ode_t* integ = context;
  real_t* xx = NV_DATA(U);
  real_t* xxd = NV_DATA(U_dot);

  // Evaluate the RHS using a solution vector with ghosts.
  memcpy(integ->U_with_ghosts, xx, sizeof(real_t) * integ->num_local_values);
  int status = integ->fe(integ->context, t, integ->U_with_ghosts, xxd);
  if (status == 0)
  {
    // Copy the local result into our solution vector.
    memcpy(xx, integ->U_with_ghosts, sizeof(real_t) * integ->num_local_values);
  }

  // Tell our observers we've computed the right hand side.
  for (int i = 0; i < integ->observers->size; ++i)
  {
    ark_ode_observer_t* obs = integ->observers->data[i];
    if (obs->fe_computed != NULL)
      obs->fe_computed(obs->context, t, xx, xxd);
  }

  STOP_FUNCTION_TIMER();
  return status;
}

// This function wraps around the user-supplied implicit RHS.
static int evaluate_fi(real_t t, N_Vector U, N_Vector U_dot, void* context)
{
  START_FUNCTION_TIMER();
  ark_ode_t* integ = context;
  real_t* xx = NV_DATA(U);
  real_t* xxd = NV_DATA(U_dot);

  // Evaluate the RHS using a solution vector with ghosts.
  memcpy(integ->U_with_ghosts, xx, sizeof(real_t) * integ->num_local_values);
  int status = integ->fi(integ->context, t, integ->U_with_ghosts, xxd);
  if (status == 0)
  {
    // Copy the local result into our solution vector.
    memcpy(xx, integ->U_with_ghosts, sizeof(real_t) * integ->num_local_values);
  }

  // Tell our observers we've computed the right hand side.
  for (int i = 0; i < integ->observers->size; ++i)
  {
    ark_ode_observer_t* obs = integ->observers->data[i];
    if (obs->fi_computed != NULL)
      obs->fi_computed(obs->context, t, xx, xxd);
  }

  STOP_FUNCTION_TIMER();
  return status;
}

static bool ark_step(void* context, real_t max_dt, real_t* t, real_t* U)
{
  START_FUNCTION_TIMER();
  ark_ode_t* integ = context;
  int status = ARK_SUCCESS;

  if (integ->first_step)
  {
    ARKodeSetInitStep(integ->arkode, max_dt);
    integ->first_step = false;
  }

  // If *t + max_dt is less than the time to which we've already integrated, 
  // we don't need to integrate; we only need to interpolate backward.
  real_t t2 = *t + max_dt;
  if (t2 > integ->t)
  {
    // Integrate to at least t -> t + max_dt.
    status = ARKode(integ->arkode, t2, integ->U, &integ->t, ARK_ONE_STEP);
    if ((status != ARK_SUCCESS) && (status != ARK_TSTOP_RETURN))
    {
      integ->status_message = get_status_message(status, integ->t);
      STOP_FUNCTION_TIMER();
      return false;
    }

    if ((t2 - *t) < (integ->t - *t))
      log_detail("ark_ode_integrator: took internal step dt = %g", integ->t - *t);
  }

  // If we integrated past t2, interpolate to t2.
  if (integ->t > t2)
  {
    status = ARKodeGetDky(integ->arkode, t2, 0, integ->U);
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
  if ((status == ARK_SUCCESS) || (status == ARK_TSTOP_RETURN))
  {
    // Copy out the solution.
    memcpy(U, NV_DATA(integ->U), sizeof(real_t) * integ->num_local_values); 
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

static bool ark_advance(void* context, real_t t1, real_t t2, real_t* U)
{
  ark_ode_t* integ = context;
  integ->t = t1;
  ARKRhsFn eval_fe = (integ->fe != NULL) ? evaluate_fe : NULL;
  ARKRhsFn eval_fi = (integ->fi != NULL) ? evaluate_fi : NULL;
  ARKodeReInit(integ->arkode, eval_fe, eval_fi, 0.0, integ->U);
  ARKodeSetStopTime(integ->arkode, t2);
  
  // Copy in the solution.
  memcpy(NV_DATA(integ->U), U, sizeof(real_t) * integ->num_local_values); 

  // Integrate.
  int status = ARKode(integ->arkode, t2, integ->U, &integ->t, ARK_NORMAL);
  
  // Clear the present status.
  if (integ->status_message != NULL)
  {
    polymec_free(integ->status_message);
    integ->status_message = NULL;
  }

  // Did it work?
  if ((status == ARK_SUCCESS) || (status == ARK_TSTOP_RETURN))
  {
    // Copy out the solution.
    memcpy(U, NV_DATA(integ->U), sizeof(real_t) * integ->num_local_values); 
    return true;
  }
  else
  {
    integ->status_message = get_status_message(status, integ->t);
    return false;
  }
}

static void ark_reset(void* context, real_t t, real_t* U)
{
  ark_ode_t* integ = context;
  integ->first_step = true;

  // Reset the preconditioner.
  if (integ->precond != NULL)
    newton_pc_reset(integ->precond, t);

  // Copy in the solution and reinitialize.
  memcpy(NV_DATA(integ->U), U, sizeof(real_t) * integ->num_local_values); 
  ARKRhsFn eval_fe = (integ->fe != NULL) ? evaluate_fe : NULL;
  ARKRhsFn eval_fi = (integ->fi != NULL) ? evaluate_fi : NULL;
  ARKodeReInit(integ->arkode, eval_fe, eval_fi, 0.0, integ->U);
  integ->t = t;
}

static void ark_dtor(void* context)
{
  ark_ode_t* integ = context;

  // Kill the preconditioner stuff.
  if (integ->precond != NULL)
    newton_pc_free(integ->precond);

  // Kill the ARKode stuff.
  polymec_free(integ->U_with_ghosts);
  N_VDestroy(integ->U);
  ARKodeFree(&integ->arkode);

  // Kill the rest.
  if (integ->status_message != NULL)
    polymec_free(integ->status_message);
  if ((integ->context != NULL) && (integ->dtor != NULL))
    integ->dtor(integ->context);
  ptr_array_free(integ->observers);
  if (integ->error_weights != NULL)
    polymec_free(integ->error_weights);
  polymec_free(integ);
}

// This function sets up the preconditioner data within the integrator.
static int set_up_preconditioner(real_t t, N_Vector U, N_Vector F,
                                 int jacobian_is_current, int* jacobian_was_updated, 
                                 real_t gamma, void* context, 
                                 N_Vector work1, N_Vector work2, N_Vector work3)
{
  ark_ode_t* integ = context;
  if (!jacobian_is_current)
  {
    // Compute the approximate Jacobian using a solution vector with ghosts.
    memcpy(integ->U_with_ghosts, NV_DATA(U), sizeof(real_t) * integ->num_local_values);
    newton_pc_setup(integ->precond, 1.0, -gamma, 0.0, t, integ->U_with_ghosts, NULL);
    *jacobian_was_updated = 1;
  }
  else
    *jacobian_was_updated = 0;
  return 0;
}

// This function solves the preconditioner equation. On input, the vector r 
// contains the right-hand side of the preconditioner system, and on output 
// it contains the solution to the system.
static int solve_preconditioner_system(real_t t, N_Vector U, N_Vector F, 
                                       N_Vector r, N_Vector z, 
                                       real_t gamma, real_t delta, 
                                       int lr, void* context, 
                                       N_Vector work)
{
  ark_ode_t* integ = context;
  
  // FIXME: Apply scaling if needed.

  // Set the preconditioner's L2 norm tolerance.
  newton_pc_set_tolerance(integ->precond, delta);

  // Solve it.
  int result = (newton_pc_solve(integ->precond, t, NV_DATA(U), NULL,
                                NV_DATA(r), NV_DATA(z))) ? 0 : 1;
  return result;
}

// Adaptor for stable_dt function
static int stable_dt(N_Vector y, real_t t, real_t* dt_stable, void* context)
{
  START_FUNCTION_TIMER();
  ark_ode_t* integ = context;
  *dt_stable = integ->stable_dt(integ->context, t, NV_DATA(y));
  STOP_FUNCTION_TIMER();
  return ARK_SUCCESS;
}

// Adaptor for J*y function
static int eval_Jy(N_Vector y, N_Vector Jy, real_t t, N_Vector U, N_Vector U_dot, void* context, N_Vector tmp)
{
  START_FUNCTION_TIMER();
  ark_ode_t* integ = context;
  real_t* xx = NV_DATA(U);
  real_t* yy = NV_DATA(y);
  real_t* my_rhs = NV_DATA(U_dot);
  real_t* temp = NV_DATA(tmp);
  real_t* jjy = NV_DATA(Jy);

  // Make sure we use ghosts.
  memcpy(integ->U_with_ghosts, xx, sizeof(real_t) * integ->num_local_values);
  int status = integ->Jy(integ->context, t, xx, my_rhs, yy, temp, jjy);

  // Tell our observers we've computed the right hand side.
  for (int i = 0; i < integ->observers->size; ++i)
  {
    ark_ode_observer_t* obs = integ->observers->data[i];
    if (obs->Jy_computed != NULL)
      obs->Jy_computed(obs->context, t, xx, my_rhs, yy, jjy);
  }

  STOP_FUNCTION_TIMER();
  return status;
}

ode_integrator_t* explicit_ark_ode_integrator_new(int order, 
                                                  MPI_Comm comm,
                                                  int num_local_values, 
                                                  int num_remote_values, 
                                                  void* context, 
                                                  int (*fe_func)(void* context, real_t t, real_t* U, real_t* fe),
                                                  real_t (*stable_dt_func)(void* context, real_t, real_t* U),
                                                  void (*dtor)(void* context))
{
  return functional_ark_ode_integrator_new(order, comm, num_local_values, 
                                           num_remote_values, context, fe_func, 
                                           NULL, stable_dt_func, dtor, 0);
}

ode_integrator_t* functional_ark_ode_integrator_new(int order, 
                                                    MPI_Comm comm,
                                                    int num_local_values, 
                                                    int num_remote_values, 
                                                    void* context, 
                                                    int (*fe_func)(void* context, real_t t, real_t* U, real_t* fe),
                                                    int (*fi_func)(void* context, real_t t, real_t* U, real_t* fi),
                                                    real_t (*stable_dt_func)(void* context, real_t, real_t* U),
                                                    void (*dtor)(void* context),
                                                    int max_anderson_accel_dim)
{
  ASSERT((max_anderson_accel_dim >= 1) || (fi_func == NULL));
  ASSERT((order >= 3) || ((order >= 2) && ((fe_func == NULL) || (fi_func == NULL))));
  ASSERT((order <= 5) || ((order <= 6) && (fi_func == NULL)));
  ASSERT(num_local_values > 0);
  ASSERT(num_remote_values >= 0);
  ASSERT((fe_func != NULL) || (fi_func != NULL));

  ark_ode_t* integ = polymec_malloc(sizeof(ark_ode_t));
  integ->comm = comm;
  integ->num_local_values = num_local_values;
  integ->num_remote_values = num_remote_values;
  integ->context = context;
  integ->fe = fe_func;
  integ->fi = fi_func;
  integ->dtor = dtor;
  integ->status_message = NULL;
  integ->max_krylov_dim = 0;
  integ->stable_dt = (fe_func != NULL) ? stable_dt_func : NULL;
  integ->Jy = NULL;
  integ->t = 0.0;
  integ->observers = ptr_array_new();
  integ->error_weights = NULL;
  integ->first_step = true;
  integ->precond = NULL;

  integ->reset_func = NULL;
  integ->setup_func = NULL;
  integ->solve_func = NULL;

  // Set up ARKode and accessories.
  integ->U = N_VNew(integ->comm, integ->num_local_values);
  integ->U_with_ghosts = polymec_malloc(sizeof(real_t) * (integ->num_local_values + integ->num_remote_values));
  integ->arkode = ARKodeCreate();
  ARKodeSetErrFile(integ->arkode, log_stream(LOG_URGENT));
  ARKodeSetOrder(integ->arkode, order);
  ARKodeSetUserData(integ->arkode, integ);
  if ((integ->fe != NULL) && (integ->stable_dt != NULL))
    ARKodeSetStabilityFn(integ->arkode, stable_dt, integ);
  if (integ->fi != NULL)
    ARKodeSetFixedPoint(integ->arkode, max_anderson_accel_dim);
  ARKRhsFn eval_fe = (integ->fe != NULL) ? evaluate_fe : NULL;
  ARKRhsFn eval_fi = (integ->fi != NULL) ? evaluate_fi : NULL;
  ARKodeInit(integ->arkode, eval_fe, eval_fi, 0.0, integ->U);
  if (fi_func == NULL)
    ARKodeSetExplicit(integ->arkode);
  else if (fe_func == NULL)
    ARKodeSetImplicit(integ->arkode);
  else
    ARKodeSetImEx(integ->arkode);

  ode_integrator_vtable vtable = {.step = ark_step, .advance = ark_advance, .reset = ark_reset, .dtor = ark_dtor};
  char name[1024];
  if (integ->fi != NULL)
    snprintf(name, 1024, "Additive Runge-Kutta (fixed-point, order %d)", order);
  else
    snprintf(name, 1024, "Explicit Runge-Kutta (order %d)", order);
  ode_integrator_t* I = ode_integrator_new(name, integ, vtable, order, 
                                           num_local_values + num_remote_values);

  // Set default tolerances.
  // relative error of 1e-4 means errors are controlled to 0.01%.
  // absolute error is set to 1 because it's completely problem dependent.
  ark_ode_integrator_set_tolerances(I, 1e-4, 1.0);

  return I;
}

ode_integrator_t* jfnk_ark_ode_integrator_new(int order,
                                              MPI_Comm comm,
                                              int num_local_values, 
                                              int num_remote_values, 
                                              void* context, 
                                              int (*fe_func)(void* context, real_t t, real_t* U, real_t* fe),
                                              int (*fi_func)(void* context, real_t t, real_t* U, real_t* fi),
                                              bool fi_is_linear,
                                              bool fi_is_time_dependent,
                                              real_t (*stable_dt_func)(void* context, real_t, real_t* U),
                                              int (*Jy_func)(void* context, real_t t, real_t* U, real_t* U_dot, real_t* y, real_t* temp, real_t* Jy),
                                              void (*dtor)(void* context),
                                              newton_pc_t* precond,
                                              jfnk_ark_krylov_t solver_type,
                                              int max_krylov_dim)
{
  ASSERT((order >= 3) || ((order >= 2) && ((fe_func == NULL) || (fi_func == NULL))));
  ASSERT((order <= 5) || ((order <= 6) && (fi_func == NULL)));
  ASSERT(num_local_values > 0);
  ASSERT(num_remote_values >= 0);
  ASSERT(fi_func != NULL);
  ASSERT(precond != NULL);
  ASSERT(max_krylov_dim > 3);
  ASSERT(!newton_pc_coefficients_fixed(precond));

  ark_ode_t* integ = polymec_malloc(sizeof(ark_ode_t));
  integ->comm = comm;
  integ->num_local_values = num_local_values;
  integ->num_remote_values = num_remote_values;
  integ->context = context;
  integ->fe = fe_func;
  integ->fi = fi_func;
  integ->dtor = dtor;
  integ->status_message = NULL;
  integ->max_krylov_dim = max_krylov_dim;
  integ->stable_dt = (fe_func != NULL) ? stable_dt_func : NULL;
  integ->Jy = Jy_func;
  integ->t = 0.0;
  integ->observers = ptr_array_new();
  integ->error_weights = NULL;
  integ->first_step = true;

  integ->reset_func = NULL;
  integ->setup_func = NULL;
  integ->solve_func = NULL;

  // Set up ARKode and accessories.
  integ->U = N_VNew(integ->comm, integ->num_local_values);
  integ->U_with_ghosts = polymec_malloc(sizeof(real_t) * (integ->num_local_values + integ->num_remote_values));
  integ->arkode = ARKodeCreate();
  ARKodeSetErrFile(integ->arkode, log_stream(LOG_URGENT));
  ARKodeSetOrder(integ->arkode, order);
  ARKodeSetUserData(integ->arkode, integ);
  if (integ->stable_dt != NULL)
    ARKodeSetStabilityFn(integ->arkode, stable_dt, integ);
  if (fi_is_linear)
    ARKodeSetLinear(integ->arkode, fi_is_time_dependent);
  ARKRhsFn eval_fe = (integ->fe != NULL) ? evaluate_fe : NULL;
  ARKRhsFn eval_fi = (integ->fi != NULL) ? evaluate_fi : NULL;
  ARKodeInit(integ->arkode, eval_fe, eval_fi, 0.0, integ->U);
  if (fe_func == NULL)
    ARKodeSetImplicit(integ->arkode);
  else
    ARKodeSetImEx(integ->arkode);

  newton_pc_side_t newton_side = newton_pc_side(precond);
  int side;
  switch (newton_side)
  {
    case NEWTON_PC_LEFT:
      side = PREC_LEFT;
      break;
    case NEWTON_PC_RIGHT:
      side = PREC_RIGHT;
      break;
    case NEWTON_PC_BOTH:
      side = PREC_BOTH;
  }

  // Set up the solver type.
  if (solver_type == JFNK_ARK_GMRES)
  {
    ARKSpgmr(integ->arkode, side, max_krylov_dim); 
    // We use modified Gram-Schmidt orthogonalization.
    ARKSpilsSetGSType(integ->arkode, MODIFIED_GS);
  }
  else if (solver_type == JFNK_ARK_FGMRES)
  {
    ARKSpfgmr(integ->arkode, side, max_krylov_dim); 
    // We use modified Gram-Schmidt orthogonalization.
    ARKSpilsSetGSType(integ->arkode, MODIFIED_GS);
  }
  else if (solver_type == JFNK_ARK_BICGSTAB)
    ARKSpbcg(integ->arkode, side, max_krylov_dim);
  else if (solver_type == JFNK_ARK_PCG)
    ARKPcg(integ->arkode, side, max_krylov_dim);
  else
    ARKSptfqmr(integ->arkode, side, max_krylov_dim);

  // Set up the Jacobian function and preconditioner.
  if (Jy_func != NULL)
    ARKSpilsSetJacTimesVecFn(integ->arkode, eval_Jy);
  integ->precond = precond;
  ARKSpilsSetPreconditioner(integ->arkode, set_up_preconditioner,
                            solve_preconditioner_system);

  ode_integrator_vtable vtable = {.step = ark_step, .advance = ark_advance, .reset = ark_reset, .dtor = ark_dtor};
  char name[1024];
  if (integ->fe != NULL)
    snprintf(name, 1024, "JFNK IMEX Additive Runge-Kutta (fixed-point, order %d)", order);
  else
    snprintf(name, 1024, "JFNK implicit Runge-Kutta (order %d)", order);
  ode_integrator_t* I = ode_integrator_new(name, integ, vtable, order,
                                           num_local_values + num_remote_values);

  // Set default tolerances.
  // relative error of 1e-4 means errors are controlled to 0.01%.
  // absolute error is set to 1 because it's completely problem dependent.
  ark_ode_integrator_set_tolerances(I, 1e-4, 1.0);

  return I;
}

//------------------------------------------------------------------------
//                    Custom ARK integrator stuff
//------------------------------------------------------------------------

static int ark_linit(ARKodeMem ark_mem)
{
  ark_ode_t* ark = ark_mem->ark_user_data;
  real_t t = ark_mem->ark_tn;
  real_t* U = NV_DATA(ark_mem->ark_ycur);
  ark->sqrtN = -1.0;
  return ark->reset_func(ark->context, t, U);
}

static int ark_lsetup(ARKodeMem ark_mem, 
                      int convfail,
                      N_Vector ypred,
                      N_Vector fpred,
                      booleantype* jcurPtr,
                      N_Vector vtemp1,
                      N_Vector vtemp2,
                      N_Vector vtemp3)
{
  ark_ode_t* ark = ark_mem->ark_user_data;
  ark_conv_status_t conv_status;
  if (convfail == ARK_NO_FAILURES)
    conv_status = ARK_CONV_NO_FAILURES;
  else if (convfail == ARK_FAIL_BAD_J)
    conv_status = ARK_CONV_BAD_J_FAILURE;
  else
    conv_status = ARK_CONV_OTHER_FAILURE;
  real_t gamma = ark_mem->ark_gamma;
  int step = (int)ark_mem->ark_nst;
  real_t t = ark_mem->ark_tn;
  bool J_updated = false;
  real_t* U_pred = NV_DATA(ypred);
  real_t* U_dot_pred = NV_DATA(fpred);
  real_t* work1 = NV_DATA(vtemp1);
  real_t* work2 = NV_DATA(vtemp2);
  real_t* work3 = NV_DATA(vtemp3);
  if (ark->sqrtN <= 0.0)
  {
    N_VConst(1.0, vtemp1);
    ark->sqrtN = sqrt(N_VDotProd(vtemp1, vtemp1));
  }
  int status = ark->setup_func(ark->context, conv_status, gamma, step, t, 
                               U_pred, U_dot_pred, &J_updated, work1, 
                               work2, work3);
  *jcurPtr = J_updated;
  return status;
}

static int ark_lsolve(ARKodeMem ark_mem, 
                      N_Vector b, 
                      N_Vector weight,
                      N_Vector ycur, 
                      N_Vector fcur)
{
  ark_ode_t* ark = ark_mem->ark_user_data;
  real_t t = ark_mem->ark_tn;
  real_t* U = NV_DATA(ycur);
  real_t* U_dot = NV_DATA(fcur);
  real_t* W = NV_DATA(weight);
  // NOTE: We multiply the residual norm tolerance by sqrt(N) to turn the 
  // NOTE: 2-norm of the residual into its WRMS norm.
  real_t res_norm_tol = 0.05 * ark->sqrtN * ark_mem->ark_eRNrm; 
  real_t* B = NV_DATA(b);
  return ark->solve_func(ark->context, t, U, U_dot, W, res_norm_tol, B);
}

static int ark_lfree(ARKodeMem ark_mem)
{
  return 0;
}

ode_integrator_t* ark_ode_integrator_new(const char* name, 
                                         int order,
                                         MPI_Comm comm,
                                         int num_local_values, 
                                         int num_remote_values, 
                                         void* context, 
                                         int (*fe_func)(void* context, real_t t, real_t* U, real_t* fe),
                                         int (*fi_func)(void* context, real_t t, real_t* U, real_t* fi),
                                         real_t (*stable_dt_func)(void* context, real_t, real_t* U),
                                         bool fi_is_linear,
                                         bool fi_is_time_dependent,
                                         int (*reset_func)(void* context, real_t t, real_t* U),
                                         int (*setup_func)(void* context, 
                                                           ark_conv_status_t conv_status, 
                                                           real_t gamma, 
                                                           int step,
                                                           real_t t, 
                                                           real_t* U_pred, 
                                                           real_t* U_dot_pred, 
                                                           bool* J_current, 
                                                           real_t* work1, real_t* work2, real_t* work3),
                                         int (*solve_func)(void* context, 
                                                           real_t t, 
                                                           real_t* U,
                                                           real_t* U_dot,
                                                           real_t* W, 
                                                           real_t res_norm_tol, 
                                                           real_t* B),
                                         void (*dtor)(void* context))
{
  ASSERT((order >= 3) || ((order >= 2) && ((fe_func == NULL) || (fi_func == NULL))));
  ASSERT((order <= 5) || ((order <= 6) && (fi_func == NULL)));
  ASSERT(num_local_values > 0);
  ASSERT(num_remote_values >= 0);
  ASSERT(fi_func != NULL);

  ark_ode_t* integ = polymec_malloc(sizeof(ark_ode_t));
  integ->comm = comm;
  integ->num_local_values = num_local_values;
  integ->num_remote_values = num_remote_values;
  integ->context = context;
  integ->fe = fe_func;
  integ->fi = fi_func;
  integ->dtor = dtor;
  integ->status_message = NULL;
  integ->max_krylov_dim = -1;
  integ->stable_dt = (fe_func != NULL) ? stable_dt_func : NULL;
  integ->Jy = NULL;
  integ->precond = NULL;
  integ->t = 0.0;
  integ->observers = ptr_array_new();
  integ->error_weights = NULL;
  integ->first_step = true;

  // Set up ARKode and accessories.
  integ->U = N_VNew(integ->comm, integ->num_local_values);
  integ->U_with_ghosts = polymec_malloc(sizeof(real_t) * (integ->num_local_values + integ->num_remote_values));
  integ->arkode = ARKodeCreate();
  ARKodeSetErrFile(integ->arkode, log_stream(LOG_URGENT));
  ARKodeSetOrder(integ->arkode, order);
  ARKodeSetUserData(integ->arkode, integ);
  if (integ->stable_dt != NULL)
    ARKodeSetStabilityFn(integ->arkode, stable_dt, integ);
  if (fi_is_linear)
    ARKodeSetLinear(integ->arkode, fi_is_time_dependent);
  ARKRhsFn eval_fe = (integ->fe != NULL) ? evaluate_fe : NULL;
  ARKRhsFn eval_fi = (integ->fi != NULL) ? evaluate_fi : NULL;
  ARKodeInit(integ->arkode, eval_fe, eval_fi, 0.0, integ->U);
  if (fe_func == NULL)
    ARKodeSetImplicit(integ->arkode);
  else
    ARKodeSetImEx(integ->arkode);

  // Set up the solver.
  integ->reset_func = reset_func;
  integ->setup_func = setup_func;
  integ->solve_func = solve_func;
  ARKodeMem ark_mem = integ->arkode;
  ark_mem->ark_linit = ark_linit;
  ark_mem->ark_lsetup = ark_lsetup;
  ark_mem->ark_lsolve = ark_lsolve;
  ark_mem->ark_lfree = ark_lfree;
  ark_mem->ark_setupNonNull = 1; // needs to be set for lsetup to be called.

  ode_integrator_vtable vtable = {.step = ark_step, 
                                  .advance = ark_advance, 
                                  .reset = ark_reset, 
                                  .dtor = ark_dtor};
  ode_integrator_t* I = ode_integrator_new(name, integ, vtable, order,
                                           num_local_values + num_remote_values);

  // Set default tolerances.
  // relative error of 1e-4 means errors are controlled to 0.01%.
  // absolute error is set to 1 because it's completely problem dependent.
  ark_ode_integrator_set_tolerances(I, 1e-4, 1.0);

  return I;
}

void* ark_ode_integrator_context(ode_integrator_t* integrator)
{
  ark_ode_t* integ = ode_integrator_context(integrator);
  return integ->context;
}

void ark_ode_integrator_set_step_controls(ode_integrator_t* integrator,
                                          real_t max_growth,
                                          real_t max_initial_growth,
                                          real_t max_convergence_cut_factor,
                                          real_t max_accuracy_cut_factor,
                                          real_t safety_factor,
                                          real_t cfl_fraction)
{
  ark_ode_t* integ = ode_integrator_context(integrator);
  ASSERT(max_growth > 1.0);
  ASSERT(max_initial_growth > 1.0);
  ASSERT(max_convergence_cut_factor > 0.0);
  ASSERT(max_convergence_cut_factor < 1.0);
  ASSERT(max_accuracy_cut_factor > 0.0);
  ASSERT(max_accuracy_cut_factor < 1.0);
  ASSERT((cfl_fraction > 0.0) || (integ->fe == NULL));
  ASSERT((cfl_fraction <= 1.0) || (integ->fe == NULL));
  ARKodeSetMaxGrowth(integ->arkode, max_growth);
  ARKodeSetMaxFirstGrowth(integ->arkode, max_initial_growth);
  ARKodeSetMaxCFailGrowth(integ->arkode, max_convergence_cut_factor);
  ARKodeSetMaxEFailGrowth(integ->arkode, max_accuracy_cut_factor);
  ARKodeSetSafetyFactor(integ->arkode, safety_factor);
  if (integ->fe != NULL)
    ARKodeSetCFLFraction(integ->arkode, cfl_fraction);
}

void ark_ode_integrator_set_predictor(ode_integrator_t* integrator, 
                                      ark_predictor_t predictor)
{
  ark_ode_t* integ = ode_integrator_context(integrator);
  int pred = 3;
  switch (predictor)
  {
    case ARK_TRIVIAL_PREDICTOR: pred = 0; break;
    case ARK_MAXORDER_PREDICTOR: pred = 1; break;
    case ARK_VARORDER_PREDICTOR: pred = 2; break;
    case ARK_CUTOFF_PREDICTOR: pred = 3; break;
    case ARK_BOOTSTRAP_PREDICTOR: pred = 4;
  }
  ARKodeSetPredictorMethod(integ->arkode, pred);
}

void ark_ode_integrator_set_max_err_test_failures(ode_integrator_t* integrator,
                                                  int max_failures)
{
  ASSERT(max_failures > 0);
  ark_ode_t* integ = ode_integrator_context(integrator);
  ARKodeSetMaxErrTestFails(integ->arkode, max_failures);
}

void ark_ode_integrator_set_max_nonlinear_iterations(ode_integrator_t* integrator,
                                                     int max_iterations)
{
  ASSERT(max_iterations > 0);
  ark_ode_t* integ = ode_integrator_context(integrator);
  ARKodeSetMaxNonlinIters(integ->arkode, max_iterations);
}

void ark_ode_integrator_set_nonlinear_convergence_coeff(ode_integrator_t* integrator,
                                                        real_t coefficient)
{
  ASSERT(coefficient > 0.0);
  ark_ode_t* integ = ode_integrator_context(integrator);
  ARKodeSetNonlinConvCoef(integ->arkode, coefficient);
}

void ark_ode_integrator_set_tolerances(ode_integrator_t* integrator,
                                      real_t relative_tol, real_t absolute_tol)
{
  ASSERT(relative_tol > 0.0);
  ASSERT(absolute_tol > 0.0);

  ark_ode_t* integ = ode_integrator_context(integrator);

  // Clear any existing error weight function.
  integ->compute_weights = NULL;
  if (integ->error_weights != NULL)
  {
    polymec_free(integ->error_weights);
    integ->error_weights = NULL;
  }

  // Set the tolerances.
  ARKodeSStolerances(integ->arkode, relative_tol, absolute_tol);
}

// Constant error weight adaptor function.
static void use_constant_weights(void* context, real_t* y, real_t* weights)
{
  ark_ode_t* integ = context;
  ASSERT(integ->error_weights != NULL);
  memcpy(weights, integ->error_weights, sizeof(real_t) * integ->num_local_values);
}

void ark_ode_integrator_set_error_weights(ode_integrator_t* integrator, real_t* weights)
{
  ark_ode_t* integ = ode_integrator_context(integrator);
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
  ark_ode_integrator_set_error_weight_function(integrator, use_constant_weights);
}

// Error weight adaptor function.
static int compute_error_weights(N_Vector y, N_Vector ewt, void* context)
{
  ark_ode_t* integ = context;
  integ->compute_weights(integ->context, NV_DATA(y), NV_DATA(ewt));

  // Check that all the weights are non-negative.
  int N = (int)NV_LOCLENGTH(y);
  for (int i = 0; i < N; ++i)
  {
    if (NV_Ith(y, i) < 0.0)
      return -1;
  }
  return 0;
}

void ark_ode_integrator_set_error_weight_function(ode_integrator_t* integrator,
                                                  void (*compute_weights)(void* context, real_t* y, real_t* weights))
{
  ark_ode_t* integ = ode_integrator_context(integrator);
  ASSERT(compute_weights != NULL);
  integ->compute_weights = compute_weights;
  ARKodeWFtolerances(integ->arkode, compute_error_weights);
}

void ark_ode_integrator_eval_fe(ode_integrator_t* integrator, real_t t, real_t* U, real_t* fe)
{
  ark_ode_t* integ = ode_integrator_context(integrator);
  if (integ->fe != NULL)
  {
    memcpy(integ->U_with_ghosts, U, sizeof(real_t) * integ->num_local_values);
    integ->fe(integ->context, t, integ->U_with_ghosts, fe);
  }
  else
    memset(fe, 0, sizeof(real_t) * integ->num_local_values);
}

void ark_ode_integrator_eval_fi(ode_integrator_t* integrator, real_t t, real_t* U, real_t* fi)
{
  ark_ode_t* integ = ode_integrator_context(integrator);
  if (integ->fi != NULL)
  {
    memcpy(integ->U_with_ghosts, U, sizeof(real_t) * integ->num_local_values);
    integ->fi(integ->context, t, integ->U_with_ghosts, fi);
  }
  else
    memset(fi, 0, sizeof(real_t) * integ->num_local_values);
}

real_t ark_ode_integrator_stable_dt(ode_integrator_t* integrator, real_t t, real_t* X)
{
  ark_ode_t* integ = ode_integrator_context(integrator);
  if (integ->stable_dt != NULL)
    return integ->stable_dt(integ->context, t, X);
  else
    return REAL_MAX;
}

newton_pc_t* ark_ode_integrator_preconditioner(ode_integrator_t* integrator)
{
  ark_ode_t* integ = ode_integrator_context(integrator);
  return integ->precond;
}

void ark_ode_integrator_get_diagnostics(ode_integrator_t* integrator, 
                                        ark_ode_integrator_diagnostics_t* diagnostics)
{
  ark_ode_t* integ = ode_integrator_context(integrator);
  diagnostics->status_message = integ->status_message; // borrowed!
  ARKodeGetIntegratorStats(integ->arkode, 
                           &diagnostics->num_steps,
                           &diagnostics->num_explicit_steps,
                           &diagnostics->num_accuracy_steps,
                           &diagnostics->num_step_attempts,
                           &diagnostics->num_fe_evaluations, 
                           &diagnostics->num_fi_evaluations, 
                           &diagnostics->num_linear_solve_setups, 
                           &diagnostics->num_error_test_failures, 
                           &diagnostics->initial_step_size, 
                           &diagnostics->last_step_size, 
                           &diagnostics->next_step_size, 
                           &diagnostics->t);
  ARKodeGetNonlinSolvStats(integ->arkode, 
                           &diagnostics->num_nonlinear_solve_iterations,
                           &diagnostics->num_nonlinear_solve_convergence_failures);
  if (integ->precond != NULL)
  {
    ARKSpilsGetNumLinIters(integ->arkode, &diagnostics->num_linear_solve_iterations);
    ARKSpilsGetNumPrecEvals(integ->arkode, &diagnostics->num_preconditioner_evaluations);
    ARKSpilsGetNumPrecSolves(integ->arkode, &diagnostics->num_preconditioner_solves);
    ARKSpilsGetNumConvFails(integ->arkode, &diagnostics->num_linear_solve_convergence_failures);
  }
  else
  {
    diagnostics->num_linear_solve_iterations = -1;
    diagnostics->num_preconditioner_evaluations = -1;
    diagnostics->num_preconditioner_solves = -1;
    diagnostics->num_linear_solve_convergence_failures = -1;
  }
}

void ark_ode_integrator_diagnostics_fprintf(ark_ode_integrator_diagnostics_t* diagnostics, 
                                            FILE* stream)
{
  if (stream == NULL) return;
  fprintf(stream, "ARK ODE integrator diagnostics:\n");
  if (diagnostics->status_message != NULL)
    fprintf(stream, "  Status: %s\n", diagnostics->status_message);
  fprintf(stream, "  Num steps: %d\n", (int)diagnostics->num_steps);
  fprintf(stream, "  Num explicit steps: %d\n", (int)diagnostics->num_explicit_steps);
  fprintf(stream, "  Num accuracy-limited steps: %d\n", (int)diagnostics->num_accuracy_steps);
  fprintf(stream, "  Num step attempts: %d\n", (int)diagnostics->num_step_attempts);
  fprintf(stream, "  Initial step size: %g\n", diagnostics->initial_step_size);
  fprintf(stream, "  Last step size: %g\n", diagnostics->last_step_size);
  fprintf(stream, "  Next step size: %g\n", diagnostics->next_step_size);
  fprintf(stream, "  Current time: %g\n", diagnostics->t);
  fprintf(stream, "  Num fe evaluations: %d\n", (int)diagnostics->num_fe_evaluations);
  fprintf(stream, "  Num fi evaluations: %d\n", (int)diagnostics->num_fi_evaluations);
  fprintf(stream, "  Num linear solve setups: %d\n", (int)diagnostics->num_linear_solve_setups);
  if (diagnostics->num_linear_solve_iterations != -1)
    fprintf(stream, "  Num linear solve iterations: %d\n", (int)diagnostics->num_linear_solve_iterations);
  if (diagnostics->num_linear_solve_convergence_failures != -1)
    fprintf(stream, "  Num linear solve convergence failures: %d\n", (int)diagnostics->num_linear_solve_convergence_failures);
  fprintf(stream, "  Num error test failures: %d\n", (int)diagnostics->num_error_test_failures);
  fprintf(stream, "  Num nonlinear solve iterations: %d\n", (int)diagnostics->num_nonlinear_solve_iterations);
  fprintf(stream, "  Num nonlinear solve convergence failures: %d\n", (int)diagnostics->num_nonlinear_solve_convergence_failures);
  if (diagnostics->num_preconditioner_evaluations != -1)
    fprintf(stream, "  Num preconditioner evaluations: %d\n", (int)diagnostics->num_preconditioner_evaluations);
  if (diagnostics->num_preconditioner_solves != -1)
    fprintf(stream, "  Num preconditioner solves: %d\n", (int)diagnostics->num_preconditioner_solves);
}

ark_ode_observer_t* ark_ode_observer_new(void* context,
                                         void (*fe_computed)(void* context, real_t t, real_t* U, real_t* U_dot),
                                         void (*fi_computed)(void* context, real_t t, real_t* U, real_t* U_dot),
                                         void (*Jy_computed)(void* context, real_t t, real_t* U, real_t* U_dot, real_t* y, real_t* Jy),
                                         void (*dtor)(void* context))
{
  ark_ode_observer_t* obs = polymec_malloc(sizeof(ark_ode_observer_t));
  obs->context = context;
  obs->fe_computed = fe_computed;
  obs->fi_computed = fi_computed;
  obs->Jy_computed = Jy_computed;
  obs->dtor = dtor;
  return obs;
}

static void ark_ode_observer_free(ark_ode_observer_t* observer)
{
  if ((observer->dtor != NULL) && (observer->context != NULL))
    observer->dtor(observer->context);
  polymec_free(observer);
}

void ark_ode_integrator_add_observer(ode_integrator_t* integrator,
                                     ark_ode_observer_t* observer)
{
  ark_ode_t* integ = ode_integrator_context(integrator);
  ptr_array_append_with_dtor(integ->observers, observer, DTOR(ark_ode_observer_free));
}

//------------------------------------------------------------------------
//                  Inexact Newton-Krylov integrator stuff
//------------------------------------------------------------------------

typedef struct
{
  MPI_Comm comm;

  // Behavior.
  void* context;
  int (*fe_func)(void* context, real_t t, real_t* U, real_t* fe);
  int (*fi_func)(void* context, real_t t, real_t* U, real_t* fi);
  real_t (*stable_dt_func)(void* context, real_t, real_t* U);
  int (*J_func)(void* context, real_t t, real_t* U, real_t* U_dot, krylov_matrix_t* J);
  void (*dtor)(void* context);

  // Linear system stuff.
  krylov_factory_t* factory;
  krylov_solver_t* solver;
  krylov_pc_t* pc;
  krylov_matrix_t* M; // Newton matrix M = (I - gamma * J).
  krylov_vector_t* X; // Solution vector.
  krylov_vector_t* B; // Right-hand side vector.
  krylov_vector_t* W; // Weights vector.

  // Metadata.
  real_t gamma_prev; // Previous value of gamma in Newton matrix.
  int prev_J_update_step; // Step during which J was previously updated.
  matrix_sparsity_t* sparsity;
  int block_size;
} ink_ark_ode_t;

static int ink_fe(void* context, real_t t, real_t* U, real_t* fe)
{
  ink_ark_ode_t* ink = context;
  return ink->fe_func(ink->context, t, U, fe);
}

static int ink_fi(void* context, real_t t, real_t* U, real_t* fi)
{
  ink_ark_ode_t* ink = context;
  return ink->fi_func(ink->context, t, U, fi);
}

static real_t ink_stable_dt(void* context, real_t t, real_t* U)
{
  ink_ark_ode_t* ink = context;
  return ink->stable_dt_func(ink->context, t, U);
}

static int ink_reset(void* context, real_t t, real_t* U)
{
  START_FUNCTION_TIMER();
  ink_ark_ode_t* ink = context;
  
  // Free any resources currently in use.
  if (ink->M != NULL)
    krylov_matrix_free(ink->M);
  if (ink->block_size > 1)
    ink->M = krylov_factory_block_matrix(ink->factory, ink->sparsity, ink->block_size);
  else
    ink->M = krylov_factory_matrix(ink->factory, ink->sparsity);
  if (ink->X != NULL)
    krylov_vector_free(ink->X);
  if (ink->B != NULL)
    krylov_vector_free(ink->B);

  // Allocate resources.
  index_t* row_dist = matrix_sparsity_row_distribution(ink->sparsity);
  ink->X = krylov_factory_vector(ink->factory, ink->comm, row_dist);
  ink->B = krylov_factory_vector(ink->factory, ink->comm, row_dist);
  ink->W = krylov_factory_vector(ink->factory, ink->comm, row_dist);

  // Set some metadata.
  ink->gamma_prev = 1.0;
  ink->prev_J_update_step = 0;

  STOP_FUNCTION_TIMER();
  return 0;
}

static int ink_setup(void* context, 
                     ark_conv_status_t conv_status, 
                     real_t gamma, 
                     int step,
                     real_t t, 
                     real_t* U_pred, 
                     real_t* U_dot_pred, 
                     bool* J_updated, 
                     real_t* work1, real_t* work2, real_t* work3)
{
  START_FUNCTION_TIMER();
  ink_ark_ode_t* ink = context;

  real_t dgamma = ABS(gamma / ink->gamma_prev - 1.0);
  log_debug("ink_ark_ode_integrator: Calculating M = I - %g * J.", gamma);
  if ((step == 0) ||
      (step > ink->prev_J_update_step + 50) ||
      ((conv_status == ARK_CONV_BAD_J_FAILURE) && (dgamma < 0.2)) || 
      (conv_status == ARK_CONV_OTHER_FAILURE))
  {
    char reason[129];
    // Call our Jacobian calculation function.
    if (step == 0)
      snprintf(reason, 128, "first step");
    else if (step > ink->prev_J_update_step + 50)
      snprintf(reason, 128, "> 50 steps since last calculation");
    else if (conv_status == ARK_CONV_BAD_J_FAILURE)
      snprintf(reason, 128, "outdated Newton matrix");
    else
      snprintf(reason, 128, "convergence failure reduced dt");
    log_debug("ink_ark_ode_integrator: Updating J (reason: %s).", reason);
    int status = ink->J_func(ink->context, t, U_pred, U_dot_pred, ink->M);
    if (status != 0)
      return status;

    // Scale the Jacobian by -gamma and add the identity matrix.
    krylov_matrix_scale(ink->M, -gamma);
    krylov_matrix_add_identity(ink->M, 1.0);

    // Use this matrix as the operator in our solver.
    krylov_solver_set_operator(ink->solver, ink->M);

    // Save some information.
    ink->prev_J_update_step = step;
    ink->gamma_prev = gamma;
    *J_updated = true;
    STOP_FUNCTION_TIMER();
    return 0;
  }
  else
  {
    if (conv_status == ARK_CONV_NO_FAILURES)
    {
      // We don't need to update J, but we still need to recompute I - gamma*J.
      krylov_matrix_add_identity(ink->M, -1.0);
      krylov_matrix_scale(ink->M, gamma/ink->gamma_prev);
      krylov_matrix_add_identity(ink->M, 1.0);
      ink->gamma_prev = gamma;
      *J_updated = false;
      STOP_FUNCTION_TIMER();
      return 0;
    }
    else
    {
      STOP_FUNCTION_TIMER();
      return 1;
    }
  }
}

static int ink_solve(void* context, 
                     real_t t, 
                     real_t* U,
                     real_t* U_dot,
                     real_t* W, 
                     real_t res_norm_tol,
                     real_t* B) 
{
  START_FUNCTION_TIMER();
  ink_ark_ode_t* ink = context;

  // Copy RHS data from B into ink->B.
  krylov_vector_copy_in(ink->B, B);

  // Copy weights into ink->W.
  krylov_vector_copy_in(ink->W, W);

  // If the 2-norm of B is less than our tolerance, return X = 0.
  real_t B_norm = krylov_vector_norm(ink->B, 2);
  if (B_norm < res_norm_tol)
  {
    log_debug("ink_ark_ode_integrator: ||B||_2 < tolerance (%g < %g), so X -> 0.", B_norm, res_norm_tol); 
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
  int num_iters;
  bool solved = krylov_solver_solve_scaled(ink->solver, ink->B, ink->W, ink->W, 
                                           ink->X, &res_norm, &num_iters);

  if (solved)
  {
    log_debug("ink_ark_ode_integrator: Solved A*X = B (||R||_2 == %g after %d iters).", res_norm, num_iters);

    // Copy solution data from ink->X into B.
    krylov_vector_copy_out(ink->X, B);
    STOP_FUNCTION_TIMER();
    return 0;
  }
  else
  {
    log_debug("ink_ark_ode_integrator: Solution to A*X = B did not converge.");
    STOP_FUNCTION_TIMER();
    return 1;
  }
}

static void ink_dtor(void* context)
{
  ink_ark_ode_t* ink = context;
  matrix_sparsity_free(ink->sparsity);
  if (ink->X != NULL)
    krylov_vector_free(ink->X);
  if (ink->B != NULL)
    krylov_vector_free(ink->B);
  if (ink->W != NULL)
    krylov_vector_free(ink->W);
  if (ink->M != NULL)
    krylov_matrix_free(ink->M);
  if (ink->pc != NULL)
    krylov_pc_free(ink->pc);
  krylov_solver_free(ink->solver);
  krylov_factory_free(ink->factory);
  if ((ink->context != NULL) && (ink->dtor != NULL))
    ink->dtor(ink->context);
}

ode_integrator_t* ink_ark_ode_integrator_new(int order, 
                                             MPI_Comm comm,
                                             krylov_factory_t* factory,
                                             matrix_sparsity_t* J_sparsity,
                                             void* context, 
                                             int (*fe_func)(void* context, real_t t, real_t* U, real_t* fe),
                                             int (*fi_func)(void* context, real_t t, real_t* U, real_t* fi),
                                             real_t (*stable_dt_func)(void* context, real_t, real_t* U),
                                             bool fi_is_linear,
                                             bool fi_is_time_dependent,
                                             int (*J_func)(void* context, real_t t, real_t* U, real_t* fi, krylov_matrix_t* J),
                                             void (*dtor)(void* context))
{
  ink_ark_ode_t* ink = polymec_malloc(sizeof(ink_ark_ode_t));
  ink->comm = comm;
  ink->context = context;
  ink->fe_func = fe_func;
  ink->fi_func = fi_func;
  ink->stable_dt_func = stable_dt_func;
  ink->J_func = J_func;
  ink->dtor = dtor;
  ink->factory = factory;
  ink->sparsity = J_sparsity;
  ink->solver = krylov_factory_gmres_solver(ink->factory, comm, 30);
  ink->pc = NULL;
  ink->M = NULL;
  ink->B = NULL;
  ink->X = NULL;
  ink->W = NULL;
  ink->block_size = 1;

  char name[1024];
  snprintf(name, 1024, "INK Additive Runge-Kutta (order %d)", order);
  int num_local_values = (int)(matrix_sparsity_num_local_rows(J_sparsity));
  ode_integrator_t* I = ark_ode_integrator_new(name, order, comm, 
                                               num_local_values, 0,
                                               ink, 
                                               (fe_func != NULL) ? ink_fe : NULL, 
                                               ink_fi, 
                                               (stable_dt_func != NULL) ? ink_stable_dt : NULL, 
                                               fi_is_linear, fi_is_time_dependent,
                                               ink_reset, ink_setup, ink_solve, ink_dtor);

  // Set default tolerances.
  // relative error of 1e-4 means errors are controlled to 0.01%.
  // absolute error is set to 1 because it's completely problem dependent.
  ark_ode_integrator_set_tolerances(I, 1e-4, 1.0);

  return I;
}

void ink_ark_ode_integrator_use_pcg(ode_integrator_t* ink_ark_ode_integ)
{
  ink_ark_ode_t* ink = ark_ode_integrator_context(ink_ark_ode_integ);
  if (ink->solver != NULL)
    krylov_solver_free(ink->solver);
  ink->solver = krylov_factory_pcg_solver(ink->factory, ink->comm);
}

void ink_ark_ode_integrator_use_gmres(ode_integrator_t* ink_ark_ode_integ,
                                      int max_krylov_dim)
{
  ink_ark_ode_t* ink = ark_ode_integrator_context(ink_ark_ode_integ);
  if (ink->solver != NULL)
    krylov_solver_free(ink->solver);
  ink->solver = krylov_factory_gmres_solver(ink->factory, ink->comm, max_krylov_dim);
}

void ink_ark_ode_integrator_use_bicgstab(ode_integrator_t* ink_ark_ode_integ)
{
  ink_ark_ode_t* ink = ark_ode_integrator_context(ink_ark_ode_integ);
  if (ink->solver != NULL)
    krylov_solver_free(ink->solver);
  ink->solver = krylov_factory_bicgstab_solver(ink->factory, ink->comm);
}

void ink_ark_ode_integrator_use_special(ode_integrator_t* ink_ark_ode_integ,
                                        const char* solver_name,
                                        string_string_unordered_map_t* options)
{
  ink_ark_ode_t* ink = ark_ode_integrator_context(ink_ark_ode_integ);
  if (ink->solver != NULL)
    krylov_solver_free(ink->solver);
  ink->solver = krylov_factory_special_solver(ink->factory, ink->comm,
                                              solver_name, options);
}

void ink_ark_ode_integrator_set_pc(ode_integrator_t* ink_ark_ode_integ,
                                   const char* pc_name, 
                                   string_string_unordered_map_t* options)
{
  ink_ark_ode_t* ink = ark_ode_integrator_context(ink_ark_ode_integ);
  if (ink->pc != NULL)
    krylov_pc_free(ink->pc);
  ink->pc = krylov_factory_preconditioner(ink->factory, ink->comm, pc_name, options);
}

void ink_ark_ode_integrator_set_block_size(ode_integrator_t* ink_ark_ode_integ,
                                           int block_size)
{
  ASSERT(block_size > 0);
  ink_ark_ode_t* ink = ark_ode_integrator_context(ink_ark_ode_integ);
  ink->block_size = block_size;
}

void* ink_ark_ode_integrator_context(ode_integrator_t* ink_ark_ode_integ)
{
  ink_ark_ode_t* ink = ark_ode_integrator_context(ink_ark_ode_integ);
  return ink->context;
}
