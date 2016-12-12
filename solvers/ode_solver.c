// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "solvers/ode_solver.h"

struct ode_solver_t 
{
  char* name;
  void* context;
  int order;
  ode_solver_vtable vtable;

  bool initialized;
  real_t current_time;

  real_t max_dt, stop_time;

  real_t* U;
  int solution_vector_size;
};

ode_solver_t* ode_solver_new(const char* name, 
                             void* context,
                             ode_solver_vtable vtable,
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

  ode_solver_t* integ = polymec_malloc(sizeof(ode_solver_t));
  integ->name = string_dup(name);
  integ->context = context;
  integ->vtable = vtable;
  integ->order = order;
  integ->solution_vector_size = solution_vector_size;
  integ->current_time = 0.0;
  integ->initialized = false;
  integ->max_dt = REAL_MAX;
  integ->stop_time = REAL_MAX;

  // We store our own copy of the solution vector if we have procedures 
  // for copying in and copying out. Otherwise, we use the actual data itself.
  if (vtable.copy_in != NULL)
    integ->U = polymec_malloc(sizeof(real_t) * integ->solution_vector_size);
  else
    integ->U = NULL;

  return integ;
}

void ode_solver_free(ode_solver_t* integ)
{
  if ((integ->context != NULL) && (integ->vtable.dtor != NULL))
    integ->vtable.dtor(integ->context);
  if (integ->U != NULL)
  {
    ASSERT(integ->vtable.copy_in != NULL);
    polymec_free(integ->U);
  }
  polymec_free(integ->name);
  polymec_free(integ);
}

char* ode_solver_name(ode_solver_t* integ)
{
  return integ->name;
}

void* ode_solver_context(ode_solver_t* integ)
{
  return integ->context;
}

int ode_solver_order(ode_solver_t* integ)
{
  return integ->order;
}

int ode_solver_solution_vector_size(ode_solver_t* integ)
{
  return integ->solution_vector_size;
}

void ode_solver_set_max_dt(ode_solver_t* integ, real_t max_dt)
{
  ASSERT(max_dt > 0.0);
  integ->max_dt = max_dt;
}

void ode_solver_set_stop_time(ode_solver_t* integ, real_t stop_time)
{
  integ->stop_time = stop_time;
}

bool ode_solver_step(ode_solver_t* integ, real_t max_dt, real_t* t, real_t* U)
{
  // Copy in data if necessary.
  if (integ->vtable.copy_in != NULL)
    integ->vtable.copy_in(integ->context, U, integ->U);
  else
    integ->U = U;

  if (!integ->initialized)
    ode_solver_reset(integ, *t, integ->U);

  // Figure out the actual maximum time.
  real_t dt = MIN(max_dt, MIN(integ->max_dt, integ->stop_time - *t));

  // Integrate.
  bool result = integ->vtable.step(integ->context, dt, t, integ->U);

  // Copy out data if necessary.
  if (integ->vtable.copy_out != NULL)
  {
    if (result)
      integ->vtable.copy_out(integ->context, integ->U, U);
  }
  else
    integ->U = NULL;

  return result;
}

bool ode_solver_advance(ode_solver_t* integ, real_t t1, real_t t2, real_t* U)
{
  // Figure out the actual end time.
  t2 = MIN(t2, integ->stop_time);

  // Copy in data if necessary.
  if (integ->vtable.copy_in != NULL)
    integ->vtable.copy_in(integ->context, U, integ->U);
  else
    integ->U = U;

  // Advance.
  bool result = integ->vtable.advance(integ->context, t1, t2, integ->U);

  // After a full integration, we must be reset.
  integ->initialized = false;

  // Copy out data if necessary.
  if (integ->vtable.copy_out != NULL)
  {
    if (result)
      integ->vtable.copy_out(integ->context, integ->U, U);
  }
  else
    integ->U = NULL;

  return result;
}

void ode_solver_reset(ode_solver_t* integ, real_t t, real_t* U)
{
  integ->current_time = t;
  if (integ->vtable.reset != NULL)
  {
    if (integ->vtable.copy_in != NULL)
    {
      if (U == integ->U)
        integ->vtable.reset(integ->context, t, U);
      else
      {
        integ->vtable.copy_in(integ->context, U, integ->U);
        integ->vtable.reset(integ->context, t, integ->U);
        integ->vtable.copy_out(integ->context, integ->U, U);
      }
    }
    else
      integ->vtable.reset(integ->context, t, U);
  }
  integ->initialized = true;
}

real_t ode_solver_current_time(ode_solver_t* integ)
{
  return integ->current_time;
}

//------------------------------------------------------------------------
//                          "Decorated" ODE solver
//------------------------------------------------------------------------

typedef struct
{
  void* context;
  ode_solver_t* base_integ;
  ode_solver_vtable vtable;
} decorated_t;

static void decorated_reset(void* context, real_t t, real_t* U)
{
  decorated_t* decorated = context;
  if (decorated->vtable.reset != NULL)
    decorated->vtable.reset(decorated->context, t, U);
  else if (decorated->base_integ->vtable.reset != NULL)
    decorated->base_integ->vtable.reset(decorated->base_integ->context, t, U);
}

static bool decorated_step(void* context, real_t max_dt, real_t* t, real_t* U)
{
  decorated_t* decorated = context;
  if (decorated->vtable.step != NULL)
    return decorated->vtable.step(decorated->context, max_dt, t, U);
  else
    return decorated->base_integ->vtable.step(decorated->base_integ->context, max_dt, t, U);
}

static bool decorated_advance(void* context, real_t t1, real_t t2, real_t* U)
{
  decorated_t* decorated = context;
  if (decorated->vtable.advance != NULL)
    return decorated->vtable.advance(decorated->context, t1, t2, U);
  else
    return decorated->base_integ->vtable.advance(decorated->base_integ->context, t1, t2, U);
}

static void decorated_copy_in(void* context, real_t* solution_data, real_t* U)
{
  decorated_t* decorated = context;
  if (decorated->vtable.copy_in != NULL)
    decorated->vtable.copy_in(decorated->context, solution_data, U);
  else if (decorated->base_integ->vtable.copy_in != NULL)
    decorated->base_integ->vtable.copy_in(decorated->base_integ->context, solution_data, U);
}

static void decorated_copy_out(void* context, real_t* U, real_t* solution_data)
{
  decorated_t* decorated = context;
  if (decorated->vtable.copy_out != NULL)
    decorated->vtable.copy_out(decorated->context, U, solution_data);
  else if (decorated->base_integ->vtable.copy_out != NULL)
    decorated->base_integ->vtable.copy_out(decorated->base_integ->context, U, solution_data);
}

static void decorated_dtor(void* context)
{
  decorated_t* decorated = context;
  ode_solver_free(decorated->base_integ);
  polymec_free(decorated);
}

ode_solver_t* decorated_ode_solver_new(ode_solver_t* base_solver,
                                       void* context,
                                       ode_solver_vtable decoration_vtable)
{
  decorated_t* decorated = polymec_malloc(sizeof(decorated_t));
  decorated->context = context;
  decorated->vtable = decoration_vtable;
  decorated->base_integ = base_solver;

  // Assemble the vtable with overrides and fallbacks.
  ode_solver_vtable vtable;
  if ((decoration_vtable.reset != NULL) || (base_solver->vtable.reset != NULL))
    vtable.reset = decorated_reset;
  vtable.step = decorated_step;
  vtable.advance = decorated_advance;
  if ((decoration_vtable.copy_in != NULL) || (base_solver->vtable.copy_in != NULL))
    vtable.copy_in = decorated_copy_in;
  if ((decoration_vtable.copy_out != NULL) || (base_solver->vtable.copy_out != NULL))
    vtable.copy_out = decorated_copy_out;
  vtable.dtor = decorated_dtor;

  int name_len = (int)strlen(base_solver->name);
  char name[name_len + 128];
  snprintf(name, name_len + 127, "Decorated %s", base_solver->name);
  ode_solver_t* I = ode_solver_new(name, decorated, vtable, 
                                   base_solver->order,
                                   base_solver->solution_vector_size);
  return I;
}

