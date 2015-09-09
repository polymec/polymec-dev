// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/unordered_map.h"
#include "multiphysics/physics_kernel.h"

struct physics_kernel_t 
{
  char* name;
  physics_kernel_type_t type;
  physics_kernel_vtable vtable;
  void* context;

  string_ptr_unordered_map_t* primary_map;
  string_ptr_unordered_map_t* secondary_map;
}; 

physics_kernel_t* physics_kernel_new(const char* name, 
                                     physics_kernel_type_t type,
                                     void* context,
                                     physics_kernel_vtable vtable)
{
  physics_kernel_t* kernel = polymec_malloc(sizeof(physics_kernel_t));
  kernel->name = string_dup(name);
  kernel->type = type;
  kernel->context = context;
  kernel->vtable = vtable;
  kernel->primary_map = string_ptr_unordered_map_new();
  kernel->secondary_map = string_ptr_unordered_map_new();
  return kernel;
}

void physics_kernel_free(physics_kernel_t* kernel)
{
  if ((kernel->vtable.dtor != NULL) && (kernel->context != NULL))
    kernel->vtable.dtor(kernel->context);
  string_ptr_unordered_map_free(kernel->primary_map);
  string_ptr_unordered_map_free(kernel->secondary_map);
  polymec_free(kernel);
}

char* physics_kernel_name(physics_kernel_t* kernel)
{
  return kernel->name;
}

physics_kernel_type_t physics_kernel_type(physics_kernel_t* kernel)
{
  return kernel->type;
}

void* physics_kernel_context(physics_kernel_t* kernel)
{
  return kernel->context;
}

int physics_kernel_solution_size(physics_kernel_t* kernel)
{
}

void physics_kernel_add_primary(physics_kernel_t* kernel,
                                const char* var_name,
                                int num_components)
{
}

bool physics_kernel_has_primary(physics_kernel_t* kernel, const char* var_name)
{
}

bool physics_kernel_next_primary(physics_kernel_t* kernel,
                                 int* pos,
                                 char** var_name,
                                 int* index,
                                 int* num_components)
{
}

void physics_kernel_add_secondary(physics_kernel_t* kernel,
                                  const char* var_name,
                                  int num_components,
                                  physics_kernel_update_function_t update)
{
}

bool physics_kernel_has_secondary(physics_kernel_t* kernel, const char* var_name)
{
}

bool physics_kernel_next_secondary(physics_kernel_t* kernel,
                                   int* pos,
                                   char** var_name,
                                   int* index,
                                   int* num_components,
                                   physics_kernel_update_function_t* update)
{
}

int physics_kernel_jacobian_block_size(physics_kernel_t* kernel)
{
}

void physics_kernel_eval_dudt(physics_kernel_t* kernel, real_t t, physics_state_t* state, real_t* dudt)
{
}

void physics_kernel_eval_residual(physics_kernel_t* kernel, real_t t, physics_state_t* state, real_t* R)
{
}

real_t physics_kernel_max_dt(physics_kernel_t* kernel, real_t t, physics_state_t* state)
{
}

void physics_kernel_compute_Jv(void* context, real_t t, physics_state_t* state, real_t* v, real_t* Jv)
{
}

void physics_kernel_compute_J(void* context, real_t t, physics_state_t* state, krylov_matrix_t* J)
{
}

