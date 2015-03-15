// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/timer.h"
#include "integrators/newton_pc.h"

struct newton_pc_t
{
  char* name;
  void* context;
  newton_pc_vtable vtable;

  // Fixed coefficients.
  bool coeffs_fixed;
  real_t alpha0, beta0, gamma0;
};

newton_pc_t* newton_pc_new(const char* name,
                           void* context,
                           newton_pc_vtable vtable)
{
  ASSERT(vtable.compute_p != NULL);
  ASSERT(vtable.solve != NULL);

  newton_pc_t* pc = polymec_malloc(sizeof(newton_pc_t));
  pc->name = string_dup(name);
  pc->context = context;
  pc->vtable = vtable;
  pc->coeffs_fixed = false;
  pc->alpha0 = pc->beta0 = pc->gamma0 = 0.0;
  
  return pc;
}

void newton_pc_free(newton_pc_t* precond)
{
  if ((precond->vtable.dtor != NULL) && (precond->context != NULL))
    precond->vtable.dtor(precond->context);
  string_free(precond->name);
  polymec_free(precond);
}

char* newton_pc_name(newton_pc_t* precond)
{
  return precond->name;
}

void* newton_pc_context(newton_pc_t* precond)
{
  return precond->context;
}

void newton_pc_setup(newton_pc_t* precond, 
                     real_t alpha, real_t beta, real_t gamma,
                     real_t t, real_t* x, real_t* xdot)
{
  START_FUNCTION_TIMER();
  log_debug("newton_pc: setting up preconditioner...");
  if (!precond->coeffs_fixed)
  {
    // Only certain combinations of alpha, beta, and gamma are allowed.
    ASSERT(((alpha == 1.0) && (beta != 0.0) && (gamma == 0.0)) || 
           ((alpha == 0.0) && (beta == 1.0)));

    precond->vtable.compute_p(precond->context, alpha, beta, gamma, t, x, xdot);
  }
  else
  {
    precond->vtable.compute_p(precond->context, precond->alpha0, precond->beta0, 
                              precond->gamma0, t, x, xdot);
  }
  STOP_FUNCTION_TIMER();
}

bool newton_pc_solve(newton_pc_t* precond, real_t* R)
{
  START_FUNCTION_TIMER();
  log_debug("newton_pc: solving preconditioner system...");
  int status = precond->vtable.solve(precond->context, R);
  STOP_FUNCTION_TIMER();
  return status;
}

void newton_pc_fix_coefficients(newton_pc_t* precond, real_t alpha0, real_t beta0, real_t gamma0)
{
  precond->coeffs_fixed = true;
  precond->alpha0 = alpha0;
  precond->beta0 = beta0;
  precond->gamma0 = gamma0;
}

void newton_pc_unfix_coefficients(newton_pc_t* precond)
{
  precond->coeffs_fixed = false;
}

bool newton_pc_coefficients_fixed(newton_pc_t* precond)
{
  return precond->coeffs_fixed;
}

