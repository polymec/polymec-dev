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

#include "integrators/ode_integrator.h"

struct ode_integrator_t 
{
  char* name;
  void* context;
  int order;
  ode_integrator_vtable vtable;

  bool initialized;
  real_t current_time;
};

ode_integrator_t* ode_integrator_new(const char* name, 
                                     void* context,
                                     ode_integrator_vtable vtable,
                                     int order)
{
  ASSERT(order > 0);
  ASSERT(vtable.step != NULL);
  ASSERT(vtable.advance != NULL);

  ode_integrator_t* integ = polymec_malloc(sizeof(ode_integrator_t));
  integ->name = string_dup(name);
  integ->context = context;
  integ->vtable = vtable;
  integ->order = order;
  integ->current_time = 0.0;
  integ->initialized = false;
  return integ;
}

void ode_integrator_free(ode_integrator_t* integ)
{
  if ((integ->context != NULL) && (integ->vtable.dtor != NULL))
    integ->vtable.dtor(integ->context);
  polymec_free(integ->name);
  polymec_free(integ);
}

char* ode_integrator_name(ode_integrator_t* integ)
{
  return integ->name;
}

void* ode_integrator_context(ode_integrator_t* integ)
{
  return integ->context;
}

int ode_integrator_order(ode_integrator_t* integ)
{
  return integ->order;
}

bool ode_integrator_step(ode_integrator_t* integ, real_t max_dt, real_t* t, real_t* x)
{
  if (!integ->initialized)
    ode_integrator_reset(integ, *t, x);

  // Integrate.
  return integ->vtable.step(integ->context, max_dt, t, x);
}

bool ode_integrator_advance(ode_integrator_t* integ, real_t t1, real_t t2, real_t* x)
{
  // Advance.
  bool result = integ->vtable.advance(integ->context, t1, t2, x);

  // After a full integration, we must be reset.
  integ->initialized = false;

  return result;
}

void ode_integrator_reset(ode_integrator_t* integrator, real_t t, real_t* x)
{
  integrator->current_time = t;
  if (integrator->vtable.reset != NULL)
    integrator->vtable.reset(integrator->context, t, x);
  integrator->initialized = true;
}

real_t ode_integrator_current_time(ode_integrator_t* integrator)
{
  return integrator->current_time;
}
