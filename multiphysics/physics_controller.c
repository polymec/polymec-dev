// Copyright (c) 2012-2016, Jeffrey N. Johnson
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
  ptr_array_append_with_dtor(controller->kernels, kernel, DTOR(physics_kernel_free));
}

physics_state_t* physics_controller_state(physics_controller_t* controller)
{
  physics_state_t* state = physics_state_new();
  
  // Add all the primary variables from the kernels.
  for (int k = 0; k < controller->kernels->size; ++k)
  {
    physics_kernel_t* kernel = controller->kernels->data[k];
    int pos = 0, index, size, num_components;
    char* var_name;
    while (physics_kernel_next_primary(kernel, &pos, &var_name, &index, &size, &num_components))
    {
      if (!physics_state_has_primary(state, var_name))
        physics_state_add_primary(state, var_name, size, num_components);
    }
  }

  // Add all the secondary variables from the kernels. Maintenance responsibility
  // is on a first-come, first-served basis.
  for (int k = 0; k < controller->kernels->size; ++k)
  {
    physics_kernel_t* kernel = controller->kernels->data[k];
    int pos = 0, index, size, num_components;
    char* var_name;
    physics_kernel_update_function_t update;
    while (physics_kernel_next_secondary(kernel, &pos, &var_name, &index, &size, &num_components, &update))
    {
      if (!physics_state_has_secondary(state, var_name))
        physics_state_add_secondary(state, var_name, size, num_components, kernel);
    }
  }

  return state;
}

void physics_controller_advance(physics_controller_t* controller,
                                real_t* t, 
                                physics_state_t* state)
{
  controller->vtable.advance(controller->context, 
                             (physics_kernel_t**)controller->kernels->data,
                             controller->kernels->size, t, state);
}

typedef struct
{
  ptr_array_t* integrators;
} op_split_ctrl_t;

static void op_split_advance(void* context, 
                             physics_kernel_t** kernels, 
                             int num_kernels,
                             real_t* t,
                             physics_state_t* state)
{
}

static void op_split_dtor(void* context)
{
  op_split_ctrl_t* c = context;
  ptr_array_free(c->integrators);
  polymec_free(c);
}

physics_controller_t* operator_split_physics_controller(physics_kernel_t** kernels,
                                                        ode_integrator_t** integrators,
                                                        int num_kernels)
{
  op_split_ctrl_t* controller = polymec_malloc(sizeof(op_split_ctrl_t));
  controller->integrators = ptr_array_new();
  for (int i = 0; i < num_kernels; ++i)
    ptr_array_append_with_dtor(controller->integrators, integrators[i], DTOR(ode_integrator_free));

  physics_controller_vtable vtable = {.advance = op_split_advance,
                                      .dtor = op_split_dtor};
  physics_controller_t* c = physics_controller_new("Operator-split", controller, vtable);
  for (int i = 0; i < num_kernels; ++i)
    physics_controller_add_kernel(c, kernels[i]);
  return c;
}

typedef struct
{
  ode_integrator_t* integrator;
} fc_ctrl_t;

static void fc_advance(void* context, 
                       physics_kernel_t** kernels, 
                       int num_kernels,
                       real_t* t,
                       physics_state_t* state)
{
}

static void fc_dtor(void* context)
{
  fc_ctrl_t* c = context;
  ode_integrator_free(c->integrator);
  polymec_free(c);
}

physics_controller_t* fully_coupled_physics_controller(physics_kernel_t** kernels,
                                                       int num_kernels,
                                                       ode_integrator_t* integrator)
{
  fc_ctrl_t* controller = polymec_malloc(sizeof(fc_ctrl_t));
  controller->integrator = integrator;
  physics_controller_vtable vtable = {.advance = fc_advance,
                                      .dtor = fc_dtor};
  physics_controller_t* c = physics_controller_new("Fully-coupled", controller, vtable);
  for (int i = 0; i < num_kernels; ++i)
    physics_controller_add_kernel(c, kernels[i]);
  return c;
}
