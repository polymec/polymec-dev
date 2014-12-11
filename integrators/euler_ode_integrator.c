// Copyright (c) 2012-2014, Jeffrey N. Johnson
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

#include "integrators/euler_ode_integrator.h"

typedef struct
{
  real_t alpha;
  MPI_Comm comm;

  // Context and v-table.
  void* context; 
  int (*rhs)(void* context, real_t t, real_t* x, real_t* xdot);
  void (*dtor)(void* context);

  // Solution data -- iterates.
  int num_local_values, num_remote_values;
  real_t *x_old, *x_new, *x_mid, *f;

  // Iteration stuff.
  int max_iters;
  real_t abs_tol, rel_tol;
} euler_ode_t;

static bool euler_step(void* context, real_t max_dt, real_t* t, real_t* x)
{
  euler_ode_t* integ = context;
  real_t alpha = integ->alpha;
  int N_local = integ->num_local_values;

  if (alpha == 0.0)
  {
    // No implicitness here! Just compute the RHS and mash it into 
    // our solution.
    int status = integ->rhs(integ->context, *t, x, integ->f);
    if (status != 0)
      return false;
    for (int i = 0; i < N_local; ++i)
      x[i] += max_dt * integ->f[i];
    *t += max_dt;
  }
  else
  {
    // Copy the solution into place.
    memcpy(integ->x_new, x, sizeof(real_t) * N_local);

    // Okay, let's do this iteration thing.
    bool converged = false;
    for (int i = 0; i < integ->max_iters; ++i)
    {
      // Stash the old iterate for x_new in x_old.
      memcpy(integ->x_old, integ->x_new, sizeof(real_t) * N_local);

      // Form x_mid from our implicitness prescription.
      for (int j = 0; j < N_local; ++j)
        integ->x_mid[j] = (1.0 - alpha) * x[j] + alpha * integ->x_new[j]; 

      // Stuff it into the RHS function.
      real_t tt = *t + integ->alpha * max_dt;
      int status = integ->rhs(integ->context, tt, integ->x_mid, integ->f);
      if (status != 0) break; // Error, not converged.

      // Update x_new: x_new <- x_orig + max_dt * RHS.
      real_t max_rel_change = 0.0, max_abs_change = 0.0;
      for (int j = 0; j < N_local; ++j)
      {
        integ->x_new[j] = x[j] + max_dt * integ->f[j];
        real_t x_new = integ->x_new[j], x_old = integ->x_old[j];
        max_abs_change = MAX(max_abs_change, fabs(x_new - x_old));
        if ((x_old == 0.0) && (x_new == 0.0))
          max_rel_change = 0.0;
        else if ((x_old == 0.0) || (x_new == 0.0))
        {
          real_t ref_value = MAX(fabs(x_old), fabs(x_new));
          max_rel_change = MAX(max_rel_change, fabs((x_new - x_old)/ref_value));
        }
        else
          max_rel_change = MAX(max_rel_change, fabs((x_new - x_old)/x_old));
      }

      // Check for convergence.
      if ((max_rel_change < integ->rel_tol) && (max_abs_change < integ->abs_tol))
      {
        memcpy(x, integ->x_new, sizeof(real_t) * N_local);
        converged = true;
        *t += max_dt;
        break;
      }
    }
    return converged;
  }
}

static bool euler_advance(void* context, real_t t1, real_t t2, real_t* x)
{
  // Well, if you're not going to be picky, we're going to be greedy. After 
  // all, the backward Euler method is L-stable, so you're probably only 
  // calling this if you're using alpha == 1. And if you're not, it's time 
  // to learn from your mistakes.
  real_t max_dt = t2 - t1; // Take one honkin' step.
  real_t t = t1;
  return euler_step(context, max_dt, &t, x);
}

static void euler_dtor(void* context)
{
  euler_ode_t* integ = context;
  polymec_free(integ->f);
  polymec_free(integ->x_new);
  polymec_free(integ->x_mid);
  polymec_free(integ->x_old);
  if ((integ->context != NULL) && (integ->dtor != NULL))
    integ->dtor(integ->context);
  polymec_free(integ);
}

ode_integrator_t* functional_euler_ode_integrator_new(real_t alpha, 
                                                      MPI_Comm comm, 
                                                      int num_local_values,
                                                      int num_remote_values,
                                                      void* context, 
                                                      int (*rhs)(void* context, real_t t, real_t* x, real_t* xdot),
                                                      void (*dtor)(void* context))
{
  ASSERT(alpha >= 0.0);
  ASSERT(alpha <= 1.0);
  ASSERT(num_local_values > 0);
  ASSERT(num_remote_values >= 0);
  ASSERT(rhs != NULL);

  euler_ode_t* integ = polymec_malloc(sizeof(euler_ode_t));
  integ->alpha = alpha;
  integ->comm = comm;
  integ->num_local_values = num_local_values;
  integ->num_remote_values = num_remote_values;
  integ->context = context;
  integ->rhs = rhs;
  integ->dtor = dtor;

  integ->x_old = polymec_malloc(sizeof(real_t) * (num_local_values + num_remote_values));
  integ->x_mid = polymec_malloc(sizeof(real_t) * (num_local_values + num_remote_values));
  integ->x_new = polymec_malloc(sizeof(real_t) * (num_local_values + num_remote_values));
  integ->f = polymec_malloc(sizeof(real_t) * num_local_values); // no ghosts here!

  ode_integrator_vtable vtable = {.step = euler_step, 
                                  .advance = euler_advance, 
                                  .dtor = euler_dtor};

  int order = 1;
  if (alpha == 0.5) 
    order = 2;

  char name[1024];
  snprintf(name, 1024, "Functional Euler integrator(alpha = %g)", alpha);
  ode_integrator_t* I = ode_integrator_new(name, integ, vtable, order);

  // Set default iteration criteria.
  euler_ode_integrator_set_max_iterations(I, 100);
  euler_ode_integrator_set_tolerances(I, 1e-4, 1.0);

  return I;
}

void euler_ode_integrator_set_max_iterations(ode_integrator_t* integrator,
                                             int max_iters)
{
  ASSERT(max_iters > 0);

  euler_ode_t* integ = ode_integrator_context(integrator);
  integ->max_iters = max_iters;
}

void euler_ode_integrator_set_tolerances(ode_integrator_t* integrator,
                                         real_t relative_tol, 
                                         real_t absolute_tol)
{
  ASSERT(relative_tol > 0.0);
  ASSERT(absolute_tol > 0.0);

  euler_ode_t* integ = ode_integrator_context(integrator);
  integ->rel_tol = relative_tol;
  integ->abs_tol = absolute_tol;
}

