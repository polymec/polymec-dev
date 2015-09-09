// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/unordered_map.h"
#include "multiphysics/physics_kernel.h"

typedef struct
{
  int index, num_components;
  physics_kernel_update_function_t update;
} var_t;

struct physics_kernel_t 
{
  // Polymorphic data.
  char* name;
  physics_kernel_type_t type;
  physics_kernel_vtable vtable;
  void* context;

  // Variable metdata.
  string_ptr_unordered_map_t* primary_map;
  string_ptr_unordered_map_t* secondary_map;

  // Working data.
  real_t *u, *w;
}; 

physics_kernel_t* physics_kernel_new(const char* name, 
                                     physics_kernel_type_t type,
                                     void* context,
                                     physics_kernel_vtable vtable)
{
  ASSERT(vtable.primary_size != NULL);
  ASSERT(vtable.secondary_size != NULL);
  ASSERT((vtable.eval_dudt != NULL) || (vtable.eval_residual != NULL));
  ASSERT((vtable.eval_dudt != NULL) || (type == PHYSICS_KERNEL_ELLIPTIC));
  ASSERT((vtable.eval_residual != NULL) || (type != PHYSICS_KERNEL_ELLIPTIC));
  ASSERT((vtable.max_dt != NULL) || (type != PHYSICS_KERNEL_HYPERBOLIC));

  physics_kernel_t* kernel = polymec_malloc(sizeof(physics_kernel_t));
  kernel->name = string_dup(name);
  kernel->type = type;
  kernel->context = context;
  kernel->vtable = vtable;
  kernel->primary_map = string_ptr_unordered_map_new();
  kernel->secondary_map = string_ptr_unordered_map_new();
  kernel->u = NULL;
  kernel->w = NULL;
  return kernel;
}

void physics_kernel_free(physics_kernel_t* kernel)
{
  if ((kernel->vtable.dtor != NULL) && (kernel->context != NULL))
    kernel->vtable.dtor(kernel->context);
  string_ptr_unordered_map_free(kernel->primary_map);
  string_ptr_unordered_map_free(kernel->secondary_map);
  if (kernel->u != NULL)
    polymec_free(kernel->u);
  if (kernel->w != NULL)
    polymec_free(kernel->w);
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

int physics_kernel_solution_vector_size(physics_kernel_t* kernel)
{
  int pos = 0;
  char* var_name;
  int index, size, num_comp;
  int sol_size = 0;
  while (physics_kernel_next_primary(kernel, &pos, &var_name, &index, &size, &num_comp))
    sol_size += size * num_comp;
  return sol_size;
}

void physics_kernel_add_primary(physics_kernel_t* kernel,
                                const char* var_name,
                                int num_components)
{
  ASSERT(num_components > 0);

  var_t* var = polymec_malloc(sizeof(var_t));
  var->index = 0; // FIXME: need a strategy for indexing.
  var->num_components = num_components;
  string_ptr_unordered_map_insert_with_kv_dtors(kernel->primary_map, (char*)var_name, var, string_free, polymec_free);
}

bool physics_kernel_has_primary(physics_kernel_t* kernel, const char* var_name)
{
  return string_ptr_unordered_map_contains(kernel->primary_map, (char*)var_name);
}

int physics_kernel_primary_size(physics_kernel_t* kernel, const char* var_name)
{
  if (physics_kernel_has_primary(kernel, var_name))
    return kernel->vtable.primary_size(kernel->context, var_name);
  else
    return -1;
}

bool physics_kernel_next_primary(physics_kernel_t* kernel,
                                 int* pos,
                                 char** var_name,
                                 int* index,
                                 int* size,
                                 int* num_components)
{
  void* entry;
  bool result = string_ptr_unordered_map_next(kernel->primary_map, pos, var_name, &entry);
  if (result)
  {
    var_t* var = entry;
    *index = var->index;
    *size = kernel->vtable.primary_size(kernel, *var_name);
    *num_components = var->num_components;
  }
  return result;
}

static int physics_kernel_secondary_vector_size(physics_kernel_t* kernel)
{
  int pos = 0;
  char* var_name;
  int index, size, num_comp;
  int sec_size = 0;
  physics_kernel_update_function_t update;
  while (physics_kernel_next_secondary(kernel, &pos, &var_name, &index, &size, &num_comp, &update))
    sec_size += size * num_comp;
  return sec_size;
}

void physics_kernel_add_secondary(physics_kernel_t* kernel,
                                  const char* var_name,
                                  int num_components,
                                  physics_kernel_update_function_t update)
{
  ASSERT(num_components > 0);
  ASSERT(update != NULL);

  var_t* var = polymec_malloc(sizeof(var_t));
  var->index = 0; // FIXME: Need a strategy for indexing.
  var->num_components = num_components;
  var->update = update;
  string_ptr_unordered_map_insert_with_kv_dtors(kernel->secondary_map, (char*)var_name, var, string_free, polymec_free);
}

bool physics_kernel_has_secondary(physics_kernel_t* kernel, const char* var_name)
{
  return string_ptr_unordered_map_contains(kernel->secondary_map, (char*)var_name);
}

int physics_kernel_secondary_size(physics_kernel_t* kernel, const char* var_name)
{
  if (physics_kernel_has_secondary(kernel, var_name))
    return kernel->vtable.secondary_size(kernel->context, var_name);
  else
    return -1;
}

bool physics_kernel_next_secondary(physics_kernel_t* kernel,
                                   int* pos,
                                   char** var_name,
                                   int* index,
                                   int* size,
                                   int* num_components,
                                   physics_kernel_update_function_t* update)
{
  void* entry;
  bool result = string_ptr_unordered_map_next(kernel->primary_map, pos, var_name, &entry);
  if (result)
  {
    var_t* var = entry;
    *index = var->index;
    *size = kernel->vtable.secondary_size(kernel, *var_name);
    *num_components = var->num_components;
    *update = var->update;
  }
  return result;
}

int physics_kernel_jacobian_block_size(physics_kernel_t* kernel)
{
  int pos = 0;
  char* var_name;
  int index, size, num_comp;
  int block_size = 0;
  while (physics_kernel_next_primary(kernel, &pos, &var_name, &index, &size, &num_comp))
    block_size += num_comp;
  return block_size;
}

static void extract_vars_from_state(physics_kernel_t* kernel, 
                                    physics_state_t* state)
{
  kernel->u = polymec_realloc(kernel->u, physics_kernel_solution_vector_size(kernel));
  kernel->w = polymec_realloc(kernel->w, physics_kernel_secondary_vector_size(kernel));
  int pos = 0;
  char* var_name;
  int index, size, num_comp;
  while (physics_kernel_next_primary(kernel, &pos, &var_name, &index, &size, &num_comp))
    physics_state_extract_primary(state, var_name, index, kernel->u);
  pos = 0;
  physics_kernel_update_function_t update;
  while (physics_kernel_next_secondary(kernel, &pos, &var_name, &index, &size, &num_comp, &update))
    physics_state_extract_secondary(state, var_name, index, kernel->w);
}

void physics_kernel_eval_dudt(physics_kernel_t* kernel, real_t t, physics_state_t* state, real_t* dudt)
{
  extract_vars_from_state(kernel, state);
  kernel->vtable.eval_dudt(kernel->context, t, kernel->u, kernel->w, dudt);
}

void physics_kernel_eval_residual(physics_kernel_t* kernel, real_t t, physics_state_t* state, real_t* R)
{
  extract_vars_from_state(kernel, state);
  kernel->vtable.eval_dudt(kernel->context, t, kernel->u, kernel->w, R);
}

real_t physics_kernel_max_dt(physics_kernel_t* kernel, real_t t, physics_state_t* state)
{
  extract_vars_from_state(kernel, state);
  return kernel->vtable.max_dt(kernel->context, t, kernel->u, kernel->w);
}

void physics_kernel_compute_Jv(physics_kernel_t* kernel, real_t t, physics_state_t* state, real_t* v, real_t* Jv)
{
  extract_vars_from_state(kernel->context, state);
  kernel->vtable.compute_Jv(kernel->context, t, kernel->u, kernel->w, v, Jv);
}

void physics_kernel_compute_J(physics_kernel_t* kernel, real_t t, physics_state_t* state, krylov_matrix_t* J)
{
  extract_vars_from_state(kernel->context, state);
  kernel->vtable.compute_J(kernel->context, t, kernel->u, kernel->w, J);
}

