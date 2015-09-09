// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/array.h"
#include "multiphysics/physics_controller.h"

struct physics_controller_t 
{
  char* name;
  void* context;
  physics_controller_vtable vtable;

  ptr_array_t* kernels;
};

physics_controller_t* physics_controller_new(const char* name, 
                                             void* context,
                                             physics_controller_vtable vtable)
{
  ASSERT(vtable.advance != NULL);
  physics_controller_t* controller = polymec_malloc(sizeof(physics_controller_t));
  controller->name = string_dup(name);
  controller->context = context;
  controller->vtable = vtable;
  controller->kernels = ptr_array_new();
  return controller;
}

void physics_controller_free(physics_controller_t* controller)
{
  if ((controller->vtable.dtor != NULL) && (controller->context != NULL))
    controller->vtable.dtor(controller->context);
  ptr_array_free(controller->kernels);
  polymec_free(controller);
}

void* physics_controller_context(physics_controller_t* controller)
{
  return controller->context;
}

void physics_controller_add_kernel(physics_controller_t* controller,
                                   physics_kernel_t* kernel)
{
  ptr_array_append(controller->kernels, kernels);
}

physics_state_t* physics_controller_state(physics_controller_t* controller)
{
  physics_state_t* state = physics_state_new();
  
  // Add all the primary variables from the kernels.
  int pos = 0, index, num_components;
  char* var_name;


  // Add all the secondary variables from the kernels. Maintenance responsibility
  // is on a first-come, first-served basis.
  return state;
}

void physics_controller_advance(physics_controller_t* controller,
                                real_t* t, 
                                physics_state_t* state)
{
  controller->vtable.advance(controller->context, controller->kernels->data,
                             controller->kernels->size, t, state);
}

