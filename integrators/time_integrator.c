// Copyright (c) 2012-2013, Jeffrey N. Johnson
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

#include "integrators/time_integrator.h"
#include "core/sundials_helpers.h"
#include "cvode/cvode_spils.h"
#include "cvode/cvode_spgmr.h"
#include "cvode/cvode_spbcgs.h"

struct time_integrator_t 
{
  char* name;
  void* context;
  time_integrator_vtable vtable;
  int order;
  time_integrator_solver_type_t solver_type;
};

time_integrator_t* time_integrator_new(const char* name, 
                                       void* context,
                                       time_integrator_vtable vtable,
                                       int order,
                                       time_integrator_solver_type_t solver_type)
{
  ASSERT(order > 0);
  time_integrator_t* integ = malloc(sizeof(time_integrator_t));
  integ->name = string_dup(name);
  integ->context = context;
  integ->vtable = vtable;
  integ->order = order;
  integ->solver_type = solver_type;
  return integ;
}

void time_integrator_free(time_integrator_t* integrator)
{
  if ((integrator->context != NULL) && (integrator->vtable.dtor != NULL))
    integrator->vtable.dtor(integrator->context);
  free(integrator->name);
  free(integrator);
}

char* time_integrator_name(time_integrator_t* integrator)
{
  return integrator->name;
}

void* time_integrator_context(time_integrator_t* integrator)
{
  return integrator->context;
}

int time_integrator_order(time_integrator_t* integrator)
{
  return integrator->order;
}

void time_integrator_step(time_integrator_t* integrator, double t1, double t2, double* X)
{
  ASSERT(t2 > t1);
  // FIXME
}

