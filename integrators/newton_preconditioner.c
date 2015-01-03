// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <gc/gc.h>
#include "slu_ddefs.h"
#include "slu_util.h"
#include "core/linear_algebra.h"
#include "core/array_utils.h"
#include "core/sundials_helpers.h"
#include "integrators/newton_preconditioner.h"

typedef struct
{
  void* context;
  newton_preconditioner_vtable vtable;
} nlpc_t;

static void nlpc_setup(void* context)
{
  polymec_error("Do not call preconditioner_setup() on a nonlinear preconditioner.\n"
                "Instead call newton_preconditioner_setup().");
}

static bool nlpc_solve(void* context, real_t* z)
{
  nlpc_t* pc = context;
  return pc->vtable.solve(pc->context, z);
}

static void nlpc_fprintf(void* context, FILE* stream)
{
  nlpc_t* pc = context;
  pc->vtable.fprintf(pc->context, stream);
}

static void nlpc_dtor(void* context)
{
  nlpc_t* pc = context;
  if ((pc->vtable.dtor != NULL) && (pc->context != NULL))
    pc->vtable.dtor(pc->context);
  polymec_free(pc);
}

preconditioner_t* newton_preconditioner_new(const char* name,
                                            void* context,
                                            newton_preconditioner_vtable vtable)
{
  // Make sure everything is there.
  ASSERT(vtable.compute_P != NULL);
  ASSERT(vtable.solve != NULL);

  nlpc_t* pc = polymec_malloc(sizeof(nlpc_t));
  pc->context = context;
  pc->vtable = vtable;

  preconditioner_vtable pc_vtable = {.setup = nlpc_setup,
                                     .solve = nlpc_solve,
                                     .fprintf = nlpc_fprintf,
                                     .dtor = nlpc_dtor};
  return preconditioner_new(name, pc, pc_vtable);
}

void newton_preconditioner_setup(preconditioner_t* precond, 
                                 real_t alpha, real_t beta, real_t gamma,
                                 real_t t, real_t* x, real_t* xdot)
{
  // Only certain combinations of alpha, beta, and gamma are allowed.
  ASSERT(((alpha == 1.0) && (beta != 0.0) && (gamma == 0.0)) || 
         ((alpha == 0.0) && (beta == 1.0)));

  nlpc_t* pc = preconditioner_context(precond);
  pc->vtable.compute_P(pc->context, alpha, beta, gamma, t, x, xdot);
}

