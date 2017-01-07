// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_ODE_SOLVER_H
#define POLYMEC_ODE_SOLVER_H

#include "core/polymec.h"

// This class provides an abstract interface for integrating systems of 
// nonlinear differential equations. These systems usually arise in the 
// context of nonlinear partial differential equations that have been 
// made semi-discrete using the method of lines.
typedef struct ode_solver_t ode_solver_t;

// This virtual table determines the implementation of the ode solver.
typedef struct
{
  // This function resets the internal state of the solver to use the 
  // given solution data U at the given time.
  void (*reset)(void* context, real_t t, real_t* U);

  // This function takes a single step in time, updating x and t and returning 
  // true on success and false on failure.
  bool (*step)(void* context, real_t max_dt, real_t* t, real_t* U);

  // This function integrates the solution x to time max_t, returning true 
  // on success and false on failure.
  bool (*advance)(void* context, real_t t1, real_t t2, real_t* U);

  // This optional function defines how non-standard solution data is copied 
  // into the solver's solution vecor x. Define this if you are 
  // implementing a custom solver that handles data in strange and 
  // wonderful layouts.
  void (*copy_in)(void* context, real_t* solution_data, real_t* U);

  // This optional function defines how non-standard solution data is copied 
  // from the solver's solution vecor U. Define this if you have defined 
  // copy_in above.
  void (*copy_out)(void* context, real_t* U, real_t* solution_data);

  // This function destroys the state (context) when the solver 
  // is destroyed.
  void (*dtor)(void* context);

} ode_solver_vtable;

// Creates an ode_solver that uses the given context and vtable to 
// implement time integration. The accuracy of the integration is determined 
// by its order.
ode_solver_t* ode_solver_new(const char* name, 
                                     void* context,
                                     ode_solver_vtable vtable,
                                     int order,
                                     int solution_vector_size);

// Creates an ode_solver that "decorates" the given solver with 
// another context and (partially-filled) virtual table, calling any 
// methods specified in the vtable with the given context, and falling 
// back to the given solver for methods not specified.
ode_solver_t* decorated_ode_solver_new(ode_solver_t* base_solver,
                                               void* context,
                                               ode_solver_vtable decoration_vtable);

// Frees an ODE solver.
void ode_solver_free(ode_solver_t* integ);

// Returns the name of the solver (internally stored).
char* ode_solver_name(ode_solver_t* integ);

// Returns the context object for this solver.
void* ode_solver_context(ode_solver_t* integ);

// Returns the order of the integration method.
int ode_solver_order(ode_solver_t* integ);

// Returns the size of the solution vector for this solver.
int ode_solver_solution_vector_size(ode_solver_t* integ);

// Sets the maximum time step size for the next integration step.
void ode_solver_set_max_dt(ode_solver_t* integ, real_t max_dt);

// Sets the time past which the solver will not step.
void ode_solver_set_stop_time(ode_solver_t* integ, real_t stop_time);

// Integrates the given solution U in place, taking a single step with a maximum 
// size of max_dt, starting at time *t and storing the new time there as well. 
// Returns true if the step succeeded, false if it failed for some reason. 
// If a step fails, both *t and U remain unchanged.
bool ode_solver_step(ode_solver_t* integ, real_t max_dt, real_t* t, real_t* U);

// Integrates the given solution U in place from time t1 to t2, taking as 
// many steps as are needed. The amount of work done by this method depends on the 
// number of steps taken, so do not call it in contexts that require a 
// predictable amount of work. Returns true if the integration succeeded, 
// false if it failed for some reason. If a step fails, U remains unchanged.
bool ode_solver_advance(ode_solver_t* integ, real_t t1, real_t t2, real_t* U);

// Resets the internal state of the solver to use the given solution data 
// U at the given time t. It is necessary to call this function when the 
// solution data has been altered since the last step.
void ode_solver_reset(ode_solver_t* integ, real_t t, real_t* U);

// Returns the current time at which the solver sits.
real_t ode_solver_current_time(ode_solver_t* integ);

#endif

