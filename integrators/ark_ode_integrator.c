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

struct ark_ode_observer_t 
{
  void* context;
  void (*fe_computed)(void* context, real_t t, real_t* x, real_t* fe);
  void (*fi_computed)(void* context, real_t t, real_t* x, real_t* fi);
  void (*Jy_computed)(void* context, real_t t, real_t* x, real_t* rhs, real_t* y, real_t* Jy);
  void (*dtor)(void* context);
};

typedef struct
{
  MPI_Comm comm;
  int num_local_values, num_remote_values;
  void* context; 
  int (*fe)(void* context, real_t t, real_t* x, real_t* fe);
  int (*fi)(void* context, real_t t, real_t* x, real_t* fi);
  real_t (*stable_dt)(void* context, real_t t, real_t* x);
  void (*dtor)(void* context);

  // ARKode data structures.
  void* arkode;
  N_Vector x; 
  real_t* x_with_ghosts;
  real_t t;
  char* status_message; // status of most recent integration.

  bool first_step;

  // JFNK stuff.
  int max_krylov_dim;
  newton_pc_t* precond;
  int (*Jy)(void* context, real_t t, real_t* x, real_t* rhstx, real_t* y, real_t* temp, real_t* Jy);

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
static int evaluate_fe(real_t t, N_Vector x, N_Vector x_dot, void* context)
{
  START_FUNCTION_TIMER();
  ark_ode_t* integ = context;
  real_t* xx = NV_DATA(x);
  real_t* xxd = NV_DATA(x_dot);

  // Evaluate the RHS using a solution vector with ghosts.
  memcpy(integ->x_with_ghosts, xx, sizeof(real_t) * integ->num_local_values);
  int status = integ->fe(integ->context, t, integ->x_with_ghosts, xxd);
  if (status == 0)
  {
    // Copy the local result into our solution vector.
    memcpy(xx, integ->x_with_ghosts, sizeof(real_t) * integ->num_local_values);
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
static int evaluate_fi(real_t t, N_Vector x, N_Vector x_dot, void* context)
{
  START_FUNCTION_TIMER();
  ark_ode_t* integ = context;
  real_t* xx = NV_DATA(x);
  real_t* xxd = NV_DATA(x_dot);

  // Evaluate the RHS using a solution vector with ghosts.
  memcpy(integ->x_with_ghosts, xx, sizeof(real_t) * integ->num_local_values);
  int status = integ->fi(integ->context, t, integ->x_with_ghosts, xxd);
  if (status == 0)
  {
    // Copy the local result into our solution vector.
    memcpy(xx, integ->x_with_ghosts, sizeof(real_t) * integ->num_local_values);
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

static bool ark_step(void* context, real_t max_dt, real_t* t, real_t* x)
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
    status = ARKode(integ->arkode, t2, integ->x, &integ->t, ARK_ONE_STEP);
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
    status = ARKodeGetDky(integ->arkode, t2, 0, integ->x);
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
    memcpy(x, NV_DATA(integ->x), sizeof(real_t) * integ->num_local_values); 
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

static bool ark_advance(void* context, real_t t1, real_t t2, real_t* x)
{
  ark_ode_t* integ = context;
  integ->t = t1;
  ARKRhsFn eval_fe = (integ->fe != NULL) ? evaluate_fe : NULL;
  ARKRhsFn eval_fi = (integ->fi != NULL) ? evaluate_fi : NULL;
  ARKodeReInit(integ->arkode, eval_fe, eval_fi, 0.0, integ->x);
  ARKodeSetStopTime(integ->arkode, t2);
  
  // Copy in the solution.
  memcpy(NV_DATA(integ->x), x, sizeof(real_t) * integ->num_local_values); 

  // Integrate.
  int status = ARKode(integ->arkode, t2, integ->x, &integ->t, ARK_NORMAL);
  
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
    memcpy(x, NV_DATA(integ->x), sizeof(real_t) * integ->num_local_values); 
    return true;
  }
  else
  {
    integ->status_message = get_status_message(status, integ->t);
    return false;
  }
}

static void ark_reset(void* context, real_t t, real_t* x)
{
  ark_ode_t* integ = context;
  integ->first_step = true;

  // Reset the preconditioner.
  if (integ->precond != NULL)
    newton_pc_reset(integ->precond, t);

  // Copy in the solution and reinitialize.
  memcpy(NV_DATA(integ->x), x, sizeof(real_t) * integ->num_local_values); 
  ARKRhsFn eval_fe = (integ->fe != NULL) ? evaluate_fe : NULL;
  ARKRhsFn eval_fi = (integ->fi != NULL) ? evaluate_fi : NULL;
  ARKodeReInit(integ->arkode, eval_fe, eval_fi, 0.0, integ->x);
  integ->t = t;
}

static void ark_dtor(void* context)
{
  ark_ode_t* integ = context;

  // Kill the preconditioner stuff.
  if (integ->precond != NULL)
    newton_pc_free(integ->precond);

  // Kill the ARKode stuff.
  polymec_free(integ->x_with_ghosts);
  N_VDestroy(integ->x);
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
static int set_up_preconditioner(real_t t, N_Vector x, N_Vector F,
                                 int jacobian_is_current, int* jacobian_was_updated, 
                                 real_t gamma, void* context, 
                                 N_Vector work1, N_Vector work2, N_Vector work3)
{
  ark_ode_t* integ = context;
  if (!jacobian_is_current)
  {
    // Compute the approximate Jacobian using a solution vector with ghosts.
    memcpy(integ->x_with_ghosts, NV_DATA(x), sizeof(real_t) * integ->num_local_values);
    newton_pc_setup(integ->precond, 1.0, -gamma, 0.0, t, integ->x_with_ghosts, NULL);
    *jacobian_was_updated = 1;
  }
  else
    *jacobian_was_updated = 0;
  return 0;
}

// This function solves the preconditioner equation. On input, the vector r 
// contains the right-hand side of the preconditioner system, and on output 
// it contains the solution to the system.
static int solve_preconditioner_system(real_t t, N_Vector x, N_Vector F, 
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
  int result = (newton_pc_solve(integ->precond, t, NV_DATA(x), NULL,
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
static int eval_Jy(N_Vector y, N_Vector Jy, real_t t, N_Vector x, N_Vector rhs, void* context, N_Vector tmp)
{
  START_FUNCTION_TIMER();
  ark_ode_t* integ = context;
  real_t* xx = NV_DATA(x);
  real_t* yy = NV_DATA(y);
  real_t* my_rhs = NV_DATA(rhs);
  real_t* temp = NV_DATA(tmp);
  real_t* jjy = NV_DATA(Jy);

  // Make sure we use ghosts.
  memcpy(integ->x_with_ghosts, xx, sizeof(real_t) * integ->num_local_values);
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
                                                  int (*fe_func)(void* context, real_t t, real_t* x, real_t* fe),
                                                  real_t (*stable_dt_func)(void* context, real_t, real_t* x),
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
                                                    int (*fe_func)(void* context, real_t t, real_t* x, real_t* fe),
                                                    int (*fi_func)(void* context, real_t t, real_t* x, real_t* fi),
                                                    real_t (*stable_dt_func)(void* context, real_t, real_t* x),
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

  // Set up ARKode and accessories.
  integ->x = N_VNew(integ->comm, integ->num_local_values);
  integ->x_with_ghosts = polymec_malloc(sizeof(real_t) * (integ->num_local_values + integ->num_remote_values));
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
  ARKodeInit(integ->arkode, eval_fe, eval_fi, 0.0, integ->x);
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
                                              int (*fe_func)(void* context, real_t t, real_t* x, real_t* fe),
                                              int (*fi_func)(void* context, real_t t, real_t* x, real_t* fi),
                                              bool fi_is_linear,
                                              bool fi_is_time_dependent,
                                              real_t (*stable_dt_func)(void* context, real_t, real_t* x),
                                              int (*Jy_func)(void* context, real_t t, real_t* x, real_t* rhs, real_t* y, real_t* temp, real_t* Jy),
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

  // Set up ARKode and accessories.
  integ->x = N_VNew(integ->comm, integ->num_local_values);
  integ->x_with_ghosts = polymec_malloc(sizeof(real_t) * (integ->num_local_values + integ->num_remote_values));
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
  ARKodeInit(integ->arkode, eval_fe, eval_fi, 0.0, integ->x);
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
  int N = NV_LOCLENGTH(y);
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

void ark_ode_integrator_eval_fe(ode_integrator_t* integrator, real_t t, real_t* X, real_t* fe)
{
  ark_ode_t* integ = ode_integrator_context(integrator);
  if (integ->fe != NULL)
  {
    memcpy(integ->x_with_ghosts, X, sizeof(real_t) * integ->num_local_values);
    integ->fe(integ->context, t, integ->x_with_ghosts, fe);
  }
  else
    memset(fe, 0, sizeof(real_t) * integ->num_local_values);
}

void ark_ode_integrator_eval_fi(ode_integrator_t* integrator, real_t t, real_t* X, real_t* fi)
{
  ark_ode_t* integ = ode_integrator_context(integrator);
  if (integ->fi != NULL)
  {
    memcpy(integ->x_with_ghosts, X, sizeof(real_t) * integ->num_local_values);
    integ->fi(integ->context, t, integ->x_with_ghosts, fi);
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
    return FLT_MAX;
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
                                         void (*fe_computed)(void* context, real_t t, real_t* x, real_t* rhs),
                                         void (*fi_computed)(void* context, real_t t, real_t* x, real_t* rhs),
                                         void (*Jy_computed)(void* context, real_t t, real_t* x, real_t* rhs, real_t* y, real_t* Jy),
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

