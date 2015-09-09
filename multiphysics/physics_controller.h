// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_PHYSICS_CONTROLLER_H
#define POLYMEC_PHYSICS_CONTROLLER_H

#include "multiphysics/physics_state.h"
#include "multiphysics/physics_kernel.h"

// This class provides an abstract interface for a controller that advances 
// a solution in time for a set of physics kernels whose actions are applied 
// to a state.
typedef struct physics_controller_t physics_controller_t;

// This virtual table determines the implementation of the kernel.
typedef struct
{
  // This function advances the given state by applying changes from the given 
  // kernels to it, also advancing the time t.
  void (*advance)(void* context, 
                  physics_kernel_t** kernels, int num_kernels,
                  real_t* t, physics_state_t* state);

  // This function destroys the context when the controller is destroyed.
  void (*dtor)(void* context);

} physics_controller_vtable;

// Creates a new physics controller that uses the given context and vtable 
// for its behavior. This controller is not associated with any physics 
// kernels.
physics_controller_t* physics_controller_new(const char* name, 
                                             void* context,
                                             physics_controller_vtable vtable);

// Frees a physics controller.
void physics_controller_free(physics_controller_t* controller);

// Returns the context object for this controller.
void* physics_controller_context(physics_controller_t* controller);

// Adds the given physics kernel to the list of kernels controlled by this
// controller. The order in which kernels are added to a controller can affect
// the indexing of the solution vector and/or the way the integration is 
// performed on the solution.
void physics_controller_add_kernel(physics_controller_t* controller,
                                   physics_kernel_t* kernel);

// Returns a newly allocated state that is capable of containing the data 
// manipulated by the given controller, and that can be passed as an argument
// to physics_controller_advance.
physics_state_t* physics_controller_state(physics_controller_t* controller);

// Advances the solution by integrating contributions from the various 
// kernel. t is a pointer to the initial time whose value will reflect the 
// final time after the advance is complete.
void physics_controller_advance(physics_controller_t* controller,
                                real_t* t, 
                                physics_state_t* state);

#endif

