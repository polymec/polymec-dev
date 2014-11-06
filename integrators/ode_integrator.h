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

#ifndef POLYMEC_ODE_INTEGRATOR_H
#define POLYMEC_ODE_INTEGRATOR_H

#include "core/polymec.h"

// This class provides an abstract interface for integrating systems of 
// nonlinear differential equations. These systems usually arise in the 
// context of nonlinear partial differential equations that have been 
// made semi-discrete using the method of lines.
typedef struct ode_integrator_t ode_integrator_t;

// This virtual table determines the implementation of the ode integrator.
typedef struct
{
  // This function resets the internal state of the integrator to use the 
  // given solution data at the given time.
  void (*reset)(void* context, real_t t, real_t* x);

  // This function takes a single step in time, updating x and t and returning 
  // true on success and false on failure.
  bool (*step)(void* context, real_t max_dt, real_t* t, real_t* x);

  // This function integrates the solution x to time max_t, returning true 
  // on success and false on failure.
  bool (*advance)(void* context, real_t t1, real_t t2, real_t* x);

  // This function destroys the state (context) when the integrator 
  // is destroyed.
  void (*dtor)(void* context);

} ode_integrator_vtable;

// Creates an ode_integrator that uses the given context and vtable to 
// implement time integration. The accuracy of the integration is determined 
// by its order.
ode_integrator_t* ode_integrator_new(const char* name, 
                                     void* context,
                                     ode_integrator_vtable vtable,
                                     int order);

// Frees an ODE integrator.
void ode_integrator_free(ode_integrator_t* integrator);

// Returns the name of the integrator (internally stored).
char* ode_integrator_name(ode_integrator_t* integrator);

// Returns the context object for this integrator.
void* ode_integrator_context(ode_integrator_t* integrator);

// Returns the order of the integration method.
int ode_integrator_order(ode_integrator_t* integrator);

// Sets the maximum time step size for the next integration step.
void ode_integrator_set_max_dt(ode_integrator_t* integ, real_t max_dt);

// Sets the time past which the integrator will not step.
void ode_integrator_set_stop_time(ode_integrator_t* integ, real_t stop_time);

// Integrates the given solution X in place, taking a single step with a maximum 
// size of max_dt, starting at time *t and storing the new time there as well. 
// Returns true if the step succeeded, false if it failed for some reason. 
// If a step fails, both *t and X remain unchanged.
bool ode_integrator_step(ode_integrator_t* integrator, real_t max_dt, real_t* t, real_t* X);

// Integrates the given solution X in place from time t1 to t2, taking as 
// many steps as are needed. The amount of work done by this method depends on the 
// number of steps taken, so do not call it in contexts that require a 
// predictable amount of work. Returns true if the integration succeeded, 
// false if it failed for some reason. If a step fails, X remains unchanged.
bool ode_integrator_advance(ode_integrator_t* integrator, real_t t1, real_t t2, real_t* X);

// Resets the internal state of the integrator to use the given solution data 
// x at the given time t. It is necessary to call this function when the 
// solution data has been altered since the last step.
void ode_integrator_reset(ode_integrator_t* integrator, real_t t, real_t* x);

// Returns the current time at which the integrator sits.
real_t ode_integrator_current_time(ode_integrator_t* integrator);

#endif

