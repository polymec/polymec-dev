// Copyright (c) 2012-2015, Jeffrey N. Johnson
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

