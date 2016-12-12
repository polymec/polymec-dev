// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/array.h"
#include "core/timer.h"
#include "solvers/euler_ode_solver.h"

static real_t relative_difference(real_t x, real_t y)
{
  real_t diff;
  if (reals_equal(x, 0.0) && reals_equal(y, 0.0))
    diff = 0.0;
  else if (reals_equal(x, 0.0) || reals_equal(y, 0.0))
  {
    real_t ref_value = MAX(ABS(x), ABS(y));
    diff = ABS((x - y)/ref_value);
  }
  else
    diff = ABS((x - y)/x);
  return diff;
}

// These method compute the L1, L2, and Linf norms of vector x, y.
static void compute_l1_norms(MPI_Comm comm, 
                             real_t* x, 
                             real_t* y,
                             int N,
                             real_t* rel_norm,
                             real_t* abs_norm)
{
  real_t abs = 0.0, rel = 0.0;
  for (int i = 0; i < N; ++i)
  {
    real_t xi = x[i], yi = y[i];
    real_t rel_diff = relative_difference(xi, yi);
    real_t abs_diff = ABS(xi-yi);
    rel += rel_diff;
    abs += abs_diff;
  }

  real_t norms[2] = {rel, abs};
  MPI_Allreduce(MPI_IN_PLACE, norms, 2, MPI_REAL_T, MPI_SUM, comm);

  *rel_norm = norms[0];
  *abs_norm = norms[1];
}

static void compute_l2_norms(MPI_Comm comm, 
                             real_t* x, 
                             real_t* y,
                             int N,
                             real_t* rel_norm,
                             real_t* abs_norm)
{
  real_t abs = 0.0, rel = 0.0;
  for (int i = 0; i < N; ++i)
  {
    real_t xi = x[i], yi = y[i];
    real_t rel_diff = relative_difference(xi, yi);
    real_t abs_diff = ABS(xi-yi);
    rel += rel_diff*rel_diff;
    abs += abs_diff*abs_diff;
  }

  real_t norms[2] = {rel, abs};
  MPI_Allreduce(MPI_IN_PLACE, norms, 2, MPI_REAL_T, MPI_SUM, comm);

  *rel_norm = sqrt(norms[0]);
  *abs_norm = sqrt(norms[1]);
}

static void compute_linf_norms(MPI_Comm comm, 
                               real_t* x, 
                               real_t* y,
                               int N,
                               real_t* rel_norm,
                               real_t* abs_norm)
{
  real_t abs = 0.0, rel = 0.0;
  for (int i = 0; i < N; ++i)
  {
    real_t xi = x[i], yi = y[i];
    real_t rel_diff = relative_difference(xi, yi);
    real_t abs_diff = ABS(xi-yi);
    rel = MAX(rel, rel_diff);
    abs = MAX(abs, abs_diff);
  }

  real_t norms[2] = {rel, abs};
  MPI_Allreduce(MPI_IN_PLACE, norms, 2, MPI_REAL_T, MPI_MAX, comm);

  *rel_norm = norms[0];
  *abs_norm = norms[1];
}

struct euler_ode_observer_t 
{
  void* context;
  void (*rhs_computed)(void* context, real_t t, real_t* x, real_t* rhs);
  void (*dtor)(void* context);
};

typedef struct
{
  MPI_Comm comm;

  // Implicitness factor.
  real_t theta;

  // Requested time step.
  real_t dt;

  // Context and v-table.
  void* context; 
  int (*rhs)(void* context, real_t t, real_t* x, real_t* xdot);
  void (*dtor)(void* context);

  // Method for computing the norm of the absolute and relative changes
  // to the solution, for determining convergence.
  void (*compute_norms)(MPI_Comm comm, real_t* x, real_t* y, int N, real_t* abs_norm, real_t* rel_norm);

  // Solution data -- iterates.
  int num_local_values, num_remote_values;
  real_t *x_old, *x_new, *f1, *f2;

  // Newton solver (if needed).
  newton_solver_t* newton;

  // Iteration stuff.
  int max_iters;
  real_t abs_tol, rel_tol;

  // Observers.
  ptr_array_t* observers;
} euler_ode_t;

static int eval_rhs(euler_ode_t* integ, real_t t, real_t* x, real_t* rhs)
{
  int status = integ->rhs(integ->context, t, x, rhs);

  // Tell our observers we've computed the right hand side.
  for (int i = 0; i < integ->observers->size; ++i)
  {
    euler_ode_observer_t* obs = integ->observers->data[i];
    if (obs->rhs_computed != NULL)
      obs->rhs_computed(obs->context, t, x, rhs);
  }

  return status;
}

static bool euler_step(void* context, real_t max_dt, real_t* t, real_t* x)
{
  START_FUNCTION_TIMER();
  euler_ode_t* integ = context;
  integ->dt = max_dt;
  real_t theta = integ->theta;
  int N_local = integ->num_local_values;

  if (reals_equal(theta, 0.0))
  {
    // No implicitness here! Just compute the RHS and mash it into 
    // our solution.

    // Copy the solution into place.
    memcpy(integ->x_new, x, sizeof(real_t) * N_local);

    int status = eval_rhs(integ, *t, integ->x_new, integ->f1);
    if (status != 0)
    {
      log_debug("euler_ode_solver: explicit call to RHS failed.");
      STOP_FUNCTION_TIMER();
      return false;
    }
    for (int i = 0; i < N_local; ++i)
      x[i] += max_dt * integ->f1[i];
    *t += max_dt;
    STOP_FUNCTION_TIMER();
    return true;
  }
  else
  {
    real_t dt = max_dt;

    // Copy the solution into place.
    memcpy(integ->x_new, x, sizeof(real_t) * N_local);

    // Compute the right-hand function at t.
    real_t t1 = *t;
    int status = eval_rhs(integ, t1, integ->x_new, integ->f1);
    if (status != 0) 
    {
      log_debug("euler_ode_solver: implicit call to RHS at t failed.");
      STOP_FUNCTION_TIMER();
      return false;
    }

    // Okay, let's do this iteration thing.
    bool converged = false;
    for (int i = 0; i < integ->max_iters; ++i)
    {
      // Copy our latest iterate into our previous one.
      memcpy(integ->x_old, integ->x_new, sizeof(real_t) * N_local);

      // Compute the right-hand function at t + dt.
      real_t t2 = t1 + dt;
      status = eval_rhs(integ, t2, integ->x_new, integ->f2);
      if (status != 0)
      {
        log_debug("euler_ode_solver: implicit call to RHS at t + dt failed.");
        break; // Error, not converged.
      }

      // Update x_new: x_new <- x_orig + dt * RHS.
      for (int j = 0; j < N_local; ++j)
        integ->x_new[j] = x[j] + dt * ((1.0 - theta) * integ->f1[j] + theta * integ->f2[j]);

      // Compute the relative and absolute norms of the change to the solution.
      real_t rel_norm = 0.0, abs_norm = 0.0;
      integ->compute_norms(integ->comm, integ->x_old, integ->x_new, integ->num_local_values,
                           &rel_norm, &abs_norm);
       
      log_debug("euler_ode_solver: iteration %d/%d", i+1, integ->max_iters);
      log_debug("  rel_norm = %g, abs_norm = %g", rel_norm, abs_norm);

      // Check for convergence.
      if ((rel_norm < integ->rel_tol) && (abs_norm < integ->abs_tol))
      {
        memcpy(x, integ->x_new, sizeof(real_t) * N_local);
        converged = true;
        *t += dt;
        break;
      }
    }

    STOP_FUNCTION_TIMER();
    return converged;
  }
}

static int evaluate_residual(void* context, real_t t, real_t* x, real_t* R)
{
  START_FUNCTION_TIMER();
  euler_ode_t* integ = context;
  int status = eval_rhs(integ, t, x, R);
  if (status == 0)
  {
    for (int i = 0; i < integ->num_local_values; ++i)
      R[i] = x[i] - integ->x_old[i] - integ->dt * R[i];
  }
  STOP_FUNCTION_TIMER();
  return status;
}

static bool newton_euler_step(void* context, real_t max_dt, real_t* t, real_t* x)
{
  START_FUNCTION_TIMER();
  euler_ode_t* integ = context;
  integ->dt = max_dt;
  int N_local = integ->num_local_values;
  int num_iters = 0;
  memcpy(integ->x_old, x, sizeof(real_t) * N_local);
  memcpy(integ->x_new, x, sizeof(real_t) * N_local);
  
  // Make sure the preconditioner computes P = I - dt * J.
  newton_pc_t* precond = newton_solver_preconditioner(integ->newton);
  newton_pc_fix_coefficients(precond, 1.0, -integ->dt, 0.0);

  // Now solve the thing, storing the solution in integ->x_new.
  bool solved = newton_solver_solve(integ->newton, *t + integ->dt, integ->x_new, &num_iters);
  if (solved)
  {
    // Increment the time and the solution.
    *t += integ->dt;
    memcpy(x, integ->x_new, sizeof(real_t) * integ->num_local_values);
  }
  STOP_FUNCTION_TIMER();
  return solved;
}

static bool euler_advance(void* context, real_t t1, real_t t2, real_t* x)
{
  // Well, if you're not going to be picky, we're going to be greedy. After 
  // all, the backward Euler method is L-stable, so you're probably only 
  // calling this if you're using theta == 1. And if you're not, it's time 
  // to learn from your mistakes.
  real_t max_dt = t2 - t1; // Take one honkin' step.
  real_t t = t1;
  return euler_step(context, max_dt, &t, x);
}

static void euler_dtor(void* context)
{
  euler_ode_t* integ = context;
  if (integ->f1 != NULL)
    polymec_free(integ->f1);
  if (integ->f2 != NULL)
    polymec_free(integ->f2);
  polymec_free(integ->x_new);
  polymec_free(integ->x_old);
  if (integ->newton != NULL)
    newton_solver_free(integ->newton);
  ptr_array_free(integ->observers);
  if ((integ->context != NULL) && (integ->dtor != NULL))
    integ->dtor(integ->context);
  polymec_free(integ);
}

ode_solver_t* functional_euler_ode_solver_new(real_t theta, 
                                              MPI_Comm comm,
                                              int num_local_values,
                                              int num_remote_values,
                                              void* context, 
                                              int (*rhs)(void* context, real_t t, real_t* x, real_t* xdot),
                                              void (*dtor)(void* context))
{
  ASSERT(theta >= 0.0);
  ASSERT(theta <= 1.0);
  ASSERT(num_local_values > 0);
  ASSERT(num_remote_values >= 0);
  ASSERT(rhs != NULL);

  euler_ode_t* integ = polymec_malloc(sizeof(euler_ode_t));
  integ->theta = theta;
  integ->comm = comm;
  integ->num_local_values = num_local_values;
  integ->num_remote_values = num_remote_values;
  integ->context = context;
  integ->rhs = rhs;
  integ->dtor = dtor;
  integ->newton = NULL;
  integ->observers = ptr_array_new();

  integ->x_old = polymec_malloc(sizeof(real_t) * (num_local_values + num_remote_values));
  integ->x_new = polymec_malloc(sizeof(real_t) * (num_local_values + num_remote_values));
  integ->f1 = polymec_malloc(sizeof(real_t) * num_local_values); // no ghosts here!
  integ->f2 = polymec_malloc(sizeof(real_t) * num_local_values); // no ghosts here!

  ode_solver_vtable vtable = {.step = euler_step, 
                                  .advance = euler_advance, 
                                  .dtor = euler_dtor};

  int order = 1;
  if (reals_equal(theta, 0.5))
    order = 2;

  char name[1024];
  snprintf(name, 1024, "Functional Euler solver(theta = %g)", theta);
  ode_solver_t* I = ode_solver_new(name, integ, vtable, order,
                                   num_local_values + num_remote_values);

  // Set default iteration criteria.
  euler_ode_solver_set_max_iterations(I, 100);
  euler_ode_solver_set_tolerances(I, 1e-4, 1.0);
  euler_ode_solver_set_lp_convergence_norm(I, 0);

  return I;
}

ode_solver_t* jfnk_euler_ode_solver_new(MPI_Comm comm,
                                        int num_local_values,
                                        int num_remote_values,
                                        void* context,
                                        int (*rhs)(void* context, real_t t, real_t* x, real_t* xdot),
                                        void (*dtor)(void* context),
                                        newton_pc_t* precond,
                                        jfnk_newton_t solver_type,
                                        int max_krylov_dim)
{
  ASSERT(num_local_values > 0);
  ASSERT(num_remote_values >= 0);
  ASSERT(max_krylov_dim > 3);
  ASSERT(rhs != NULL);

  euler_ode_t* integ = polymec_malloc(sizeof(euler_ode_t));
  integ->theta = -REAL_MAX;
  integ->comm = comm;
  integ->num_local_values = num_local_values;
  integ->num_remote_values = num_remote_values;
  integ->context = context;
  integ->rhs = rhs;
  integ->dtor = dtor;

  // Set up the Newton solver.
  integ->newton = jfnk_newton_solver_new(comm, num_local_values, num_remote_values,
                                         integ, evaluate_residual, NULL, NULL, 
                                         NEWTON_FULL_STEP, precond, 
                                         solver_type, max_krylov_dim, 5);

  integ->observers = ptr_array_new();

  integ->x_old = polymec_malloc(sizeof(real_t) * (num_local_values + num_remote_values));
  integ->x_new = polymec_malloc(sizeof(real_t) * (num_local_values + num_remote_values));
  integ->f1 = integ->f2 = NULL;

  ode_solver_vtable vtable = {.step = newton_euler_step, 
                              .advance = euler_advance, 
                              .dtor = euler_dtor};

  int order = 1;
  ode_solver_t* I = ode_solver_new("Backward Euler Newton solver", integ, vtable, order,
                                   num_local_values + num_remote_values);

  // Set default iteration criteria.
  euler_ode_solver_set_max_iterations(I, 100);

  return I;
}

void euler_ode_solver_set_max_iterations(ode_solver_t* solver,
                                         int max_iters)
{
  ASSERT(max_iters > 0);

  euler_ode_t* integ = ode_solver_context(solver);
  if (integ->newton != NULL)
    newton_solver_set_max_iterations(integ->newton, max_iters);
  else
    integ->max_iters = max_iters;
}

void euler_ode_solver_set_tolerances(ode_solver_t* solver,
                                     real_t relative_tol, 
                                     real_t absolute_tol)
{
  ASSERT(relative_tol > 0.0);
  ASSERT(absolute_tol > 0.0);

  euler_ode_t* integ = ode_solver_context(solver);

  if (integ->newton != NULL)
    newton_solver_set_tolerances(integ->newton, absolute_tol, relative_tol);
  else
  {
    integ->rel_tol = relative_tol;
    integ->abs_tol = absolute_tol;
  }
}

void euler_ode_solver_set_lp_convergence_norm(ode_solver_t* solver,
                                              int p)
{
  ASSERT((p == 0) || (p == 1) || (p == 2));

  euler_ode_t* integ = ode_solver_context(solver);
  if (p == 0)
    integ->compute_norms = compute_linf_norms;
  else if (p == 1)
    integ->compute_norms = compute_l1_norms;
  else
    integ->compute_norms = compute_l2_norms;
}

void euler_ode_solver_set_custom_convergence_norm(ode_solver_t* solver,
                                                  void (*compute_norms)(MPI_Comm comm, real_t* x, real_t* y, int N, real_t* abs_norm, real_t* rel_norm))
{
  ASSERT(compute_norms != NULL);
  euler_ode_t* integ = ode_solver_context(solver);
  integ->compute_norms = compute_norms;
}

newton_solver_t* newton_euler_ode_solver_solver(ode_solver_t* solver)
{
  euler_ode_t* integ = ode_solver_context(solver);
  ASSERT(integ->newton != NULL);
  return integ->newton;
}
