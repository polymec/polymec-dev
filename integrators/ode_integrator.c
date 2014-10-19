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

  // CVODE data structures.
  real_t current_time, max_dt, stop_time;
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
  integ->max_dt = FLT_MAX;
  integ->stop_time = FLT_MAX;
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

void ode_integrator_set_max_dt(ode_integrator_t* integ, real_t max_dt)
{
  ASSERT(max_dt > 0);
  integ->max_dt = max_dt;
}

void ode_integrator_set_stop_time(ode_integrator_t* integ, real_t stop_time)
{
  integ->stop_time = stop_time;
}

bool ode_integrator_step(ode_integrator_t* integ, real_t* t, real_t* X)
{
  // Integrate.
  return integ->vtable.step(integ->context, integ->max_dt, integ->stop_time, t, X);
}

bool ode_integrator_advance(ode_integrator_t* integ, real_t t, real_t* X)
{
  // Advance.
  return integ->vtable.advance(integ->context, integ->max_dt, t, X);
}

