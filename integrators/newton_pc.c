// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "integrators/newton_pc.h"

struct newton_pc_t
{
  char* name;
  void* context;
  newton_pc_vtable vtable;
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
  // Only certain combinations of alpha, beta, and gamma are allowed.
  ASSERT(((alpha == 1.0) && (beta != 0.0) && (gamma == 0.0)) || 
         ((alpha == 0.0) && (beta == 1.0)));

  precond->vtable.compute_p(precond->context, alpha, beta, gamma, t, x, xdot);
}

bool newton_pc_solve(newton_pc_t* precond, real_t* R)
{
  return precond->vtable.solve(precond->context, R);
}

