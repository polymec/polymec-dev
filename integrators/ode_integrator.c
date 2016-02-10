// Copyright (c) 2012-2016, Jeffrey N. Johnson
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

  real_t max_dt, stop_time;

  real_t* x;
  int solution_vector_size;
};

ode_integrator_t* ode_integrator_new(const char* name, 
                                     void* context,
                                     ode_integrator_vtable vtable,
                                     int order,
                                     int solution_vector_size)
{
  ASSERT(order > 0);
  ASSERT(vtable.step != NULL);
  ASSERT(vtable.advance != NULL);
  ASSERT(solution_vector_size > 0);

  // Both or neithor of copy_in/copy_out must be given.
  ASSERT(((vtable.copy_in != NULL) && (vtable.copy_out != NULL)) || 
         ((vtable.copy_in == NULL) && (vtable.copy_out == NULL)));

  ode_integrator_t* integ = polymec_malloc(sizeof(ode_integrator_t));
  integ->name = string_dup(name);
  integ->context = context;
  integ->vtable = vtable;
  integ->order = order;
  integ->solution_vector_size = solution_vector_size;
  integ->current_time = 0.0;
  integ->initialized = false;
  integ->max_dt = FLT_MAX;
  integ->stop_time = FLT_MAX;

  // We store our own copy of the solution vector if we have procedures 
  // for copying in and copying out. Otherwise, we use the actual data itself.
  if (vtable.copy_in != NULL)
    integ->x = polymec_malloc(sizeof(real_t) * integ->solution_vector_size);
  else
    integ->x = NULL;

  return integ;
}

void ode_integrator_free(ode_integrator_t* integ)
{
  if ((integ->context != NULL) && (integ->vtable.dtor != NULL))
    integ->vtable.dtor(integ->context);
  if (integ->x != NULL)
  {
    ASSERT(integ->vtable.copy_in != NULL);
    polymec_free(integ->x);
  }
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

int ode_integrator_solution_vector_size(ode_integrator_t* integ)
{
  return integ->solution_vector_size;
}

void ode_integrator_set_max_dt(ode_integrator_t* integ, real_t max_dt)
{
  ASSERT(max_dt > 0.0);
  integ->max_dt = max_dt;
}

void ode_integrator_set_stop_time(ode_integrator_t* integ, real_t stop_time)
{
  integ->stop_time = stop_time;
}

bool ode_integrator_step(ode_integrator_t* integ, real_t max_dt, real_t* t, real_t* x)
{
  // Copy in data if necessary.
  if (integ->vtable.copy_in != NULL)
    integ->vtable.copy_in(integ->context, x, integ->x);
  else
    integ->x = x;

  if (!integ->initialized)
    ode_integrator_reset(integ, *t, integ->x);

  // Figure out the actual maximum time.
  real_t dt = MIN(max_dt, MIN(integ->max_dt, integ->stop_time - *t));

  // Integrate.
  bool result = integ->vtable.step(integ->context, dt, t, integ->x);

  // Copy out data if necessary.
  if (integ->vtable.copy_out != NULL)
  {
    if (result)
      integ->vtable.copy_out(integ->context, integ->x, x);
  }
  else
    integ->x = NULL;

  return result;
}

bool ode_integrator_advance(ode_integrator_t* integ, real_t t1, real_t t2, real_t* x)
{
  // Figure out the actual end time.
  t2 = MIN(t2, integ->stop_time);

  // Copy in data if necessary.
  if (integ->vtable.copy_in != NULL)
    integ->vtable.copy_in(integ->context, x, integ->x);
  else
    integ->x = x;

  // Advance.
  bool result = integ->vtable.advance(integ->context, t1, t2, integ->x);

  // After a full integration, we must be reset.
  integ->initialized = false;

  // Copy out data if necessary.
  if (integ->vtable.copy_out != NULL)
  {
    if (result)
      integ->vtable.copy_out(integ->context, integ->x, x);
  }
  else
    integ->x = NULL;

  return result;
}

void ode_integrator_reset(ode_integrator_t* integ, real_t t, real_t* x)
{
  integ->current_time = t;
  if (integ->vtable.reset != NULL)
  {
    if (integ->vtable.copy_in != NULL)
    {
      if (x == integ->x)
        integ->vtable.reset(integ->context, t, x);
      else
      {
        integ->vtable.copy_in(integ->context, x, integ->x);
        integ->vtable.reset(integ->context, t, integ->x);
        integ->vtable.copy_out(integ->context, integ->x, x);
      }
    }
    else
      integ->vtable.reset(integ->context, t, x);
  }
  integ->initialized = true;
}

real_t ode_integrator_current_time(ode_integrator_t* integ)
{
  return integ->current_time;
}

//------------------------------------------------------------------------
//                          "Decorated" ODE integrator
//------------------------------------------------------------------------

typedef struct
{
  void* context;
  ode_integrator_t* base_integ;
} decorated_t;

static void base_reset(void* context, real_t t, real_t* x)
{
  decorated_t* decorated = context;
  ode_integrator_reset(decorated->base_integ, t, x);
}
static bool base_step(void* context, real_t max_dt, real_t* t, real_t* x)
{
  decorated_t* decorated = context;
  return ode_integrator_step(decorated->base_integ, max_dt, t, x);
}

static bool base_advance(void* context, real_t t1, real_t t2, real_t* x)
{
  decorated_t* decorated = context;
  return ode_integrator_advance(decorated->base_integ, t1, t2, x);
}

static void base_copy_in(void* context, real_t* solution_data, real_t* x)
{
  decorated_t* decorated = context;
  decorated->base_integ->vtable.copy_in(decorated->base_integ->context, solution_data, x);
}

static void base_copy_out(void* context, real_t* x, real_t* solution_data)
{
  decorated_t* decorated = context;
  decorated->base_integ->vtable.copy_out(decorated->base_integ->context, x, solution_data);
}

static void decorated_dtor(void* context)
{
  decorated_t* decorated = context;
  ode_integrator_free(decorated->base_integ);
  polymec_free(decorated);
}

ode_integrator_t* decorated_ode_integrator_new(ode_integrator_t* base_integrator,
                                               void* context,
                                               ode_integrator_vtable decoration_vtable)
{
  decorated_t* decorated = polymec_malloc(sizeof(decorated_t));
  decorated->context = context;
  decorated->base_integ = base_integrator;

  // Assemble the vtable with overrides and fallbacks.
  ode_integrator_vtable vtable;
  if (decoration_vtable.reset != NULL)
    vtable.reset = decoration_vtable.reset;
  else
    vtable.reset = base_reset;
  if (decoration_vtable.step != NULL)
    vtable.step = decoration_vtable.step;
  else
    vtable.step = base_step;
  if (decoration_vtable.advance != NULL)
    vtable.advance = decoration_vtable.advance;
  else
    vtable.advance = base_advance;
  if (decoration_vtable.copy_in != NULL)
    vtable.copy_in = decoration_vtable.copy_in;
  else
    vtable.copy_in = base_copy_in;
  if (decoration_vtable.copy_out != NULL)
    vtable.copy_out = decoration_vtable.copy_out;
  else
    vtable.copy_out = base_copy_out;
  vtable.dtor = decorated_dtor;

  int name_len = strlen(base_integrator->name);
  char name[name_len + 128];
  snprintf(name, name_len + 127, "Decorated %s", base_integrator->name);
  ode_integrator_t* I = ode_integrator_new(name, decorated, vtable, 
                                           base_integrator->order,
                                           base_integrator->solution_vector_size);
  return I;
}

