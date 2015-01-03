// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

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
