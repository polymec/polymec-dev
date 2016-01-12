// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/sundials_helpers.h" // For UNIT_ROUNDOFF
#include "integrators/ssor_newton_pc.h"

typedef struct
{
  void* context;
  real_t (*f)(void* context, int i, real_t t, real_t* x);
  real_t (*DJ)(void* context, int i, real_t t, real_t* x);
  void (*dtor)(void* context);
  int num_local_values;
  real_t omega, d;
  real_t* x_pert;
} ssor_pc_t;

// This helper updates z[i] (in z) and returns the amount by which the 
// convergence norm should be incremented.
static void ssor_compute_relaxed_value(ssor_pc_t* ssor, int i, 
                                       real_t t, real_t* x, real_t* r, 
                                       real_t* z)
{
  // Compute the perturbed value of x.
  for (int i = 0; i < ssor->num_local_values; ++i)
    ssor->x_pert[i] = x[i] + ssor->d * z[i];

  // Compute a difference quotient to estimate J*z - r at i.
  real_t f_hi = ssor->f(ssor->context, i, t, ssor->x_pert);
  real_t f_lo = ssor->f(ssor->context, i, t, x);
  real_t Jz = (f_hi - f_lo) / ssor->d;
  real_t Fi = Jz - r[i];

  // Compute the diagonal element of J at i using the perturbed value of x.
  real_t DJi = ssor->DJ(ssor->context, i, t, ssor->x_pert);

  // Update zi.
  z[i] -= ssor->omega * Fi / DJi;
}

static bool ssor_newton_pc_solve(void* context, 
                                 real_t t, real_t* x, real_t* xdot, real_t tolerance,
                                 real_t* r, real_t* z, real_t* error_L2_norm)
{
  ssor_pc_t* ssor = context;
  int N = ssor->num_local_values;

  // Initialize z to zero.
  memset(z, 0, sizeof(real_t) * N);

  // Traverse the indices upward.
  for (int i = 0; i < N; ++i)
    ssor_compute_relaxed_value(ssor, i, t, x, r, z);

  // Traverse the indices downward.
  for (int i = N-1; i >= 0; --i)
    ssor_compute_relaxed_value(ssor, i, t, x, r, z);

  // Compute the error L2 norm.
  *error_L2_norm = 0.0;
  for (int i = 0; i < N; ++i)
  {
    real_t f_hi = ssor->f(ssor->context, i, t, ssor->x_pert);
    real_t f_lo = ssor->f(ssor->context, i, t, x);
    real_t Jz = (f_hi - f_lo) / ssor->d;
    real_t Fi = Jz - r[i];
    *error_L2_norm += Fi*Fi;
  }
  *error_L2_norm = sqrt(*error_L2_norm);

  return (*error_L2_norm < tolerance);
}

static void ssor_newton_pc_free(void* context)
{
  ssor_pc_t* ssor = context;
  polymec_free(ssor->x_pert);
  polymec_free(ssor);
}

newton_pc_t* ssor_newton_pc_new(void* context,
                                real_t (*f)(void* context, int i, real_t t, real_t* x),
                                real_t (*DJ)(void* context, int i, real_t t, real_t* x),
                                void (*dtor)(void* context),
                                newton_pc_side_t side,
                                int num_local_values,
                                real_t omega)
{
  ASSERT(omega > 0.0);
  ASSERT(omega < 2.0);
  ASSERT(f != NULL);
  ASSERT(DJ != NULL);
  ASSERT(num_local_values > 0);

  ssor_pc_t* ssor = polymec_malloc(sizeof(ssor_pc_t));
  ssor->context = context;
  ssor->f = f;
  ssor->DJ = DJ;
  ssor->dtor = dtor;
  ssor->num_local_values = num_local_values;
  ssor->x_pert = polymec_malloc(sizeof(real_t) * num_local_values);
  ssor->omega = omega;

  // Machine roundoff(ish).
  ssor->d = sqrt(UNIT_ROUNDOFF);

  newton_pc_vtable vtable = {.solve = ssor_newton_pc_solve,
                             .dtor = ssor_newton_pc_free};
  return newton_pc_new("SSOR preconditioner", ssor, vtable, side);
}
                                        
bool newton_pc_is_ssor_newton_pc(newton_pc_t* newton_pc)
{
  return (strcmp(newton_pc_name(newton_pc), "SSOR preconditioner") == 0);
}

void ssor_newton_pc_set_relaxation_factor(newton_pc_t* ssor_newton_pc,
                                          real_t omega)
{
  ASSERT(newton_pc_is_ssor_newton_pc(ssor_newton_pc));
  ssor_pc_t* ssor = newton_pc_context(ssor_newton_pc);
  ssor->omega = omega;
}

