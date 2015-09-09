// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/array.h"
#include "core/string_utils.h"
#include "core/unordered_map.h"
#include "multiphysics/physics_state.h"
#include "multiphysics/physics_kernel.h"

typedef struct
{
  int index, size, num_components;
} primary_var_t;

typedef struct
{
  int size, num_components;
  real_t* data;
  physics_kernel_t* kernel;
  string_array_t* deps;
} secondary_var_t;

static void secondary_var_free(void* data)
{
  secondary_var_t* var = data;
  polymec_free(var->data);
  string_array_free(var->deps);
  polymec_free(var);
}

struct physics_state_t 
{
  // Information about the solution vector for the total system.
  real_t* solution; 
  int solution_size; 
  int max_index;

  // Primary and secondary variable registries.
  string_ptr_unordered_map_t* primary_map;   // <-- only metadata is stored here
  string_ptr_unordered_map_t* secondary_map; // <-- metadata + data is stored here
};

physics_state_t* physics_state_new()
{
  physics_state_t* state = polymec_malloc(sizeof(physics_state_t));
  state->solution = NULL;
  state->solution_size = 0;
  state->max_index = 0;
  state->primary_map = string_ptr_unordered_map_new();
  state->secondary_map = string_ptr_unordered_map_new();
  return state;
}

void physics_state_free(physics_state_t* state)
{
  if (state->solution != NULL)
    polymec_free(state->solution);
  string_ptr_unordered_map_free(state->primary_map);
  string_ptr_unordered_map_free(state->secondary_map);
  polymec_free(state);
}

static void physics_state_add_primary_with_index(physics_state_t* state,
                                                 const char* var_name,
                                                 int index,
                                                 int size,
                                                 int num_components)
{
  ASSERT(index >= 0);
  ASSERT(size > 0);
  ASSERT(num_components > 0);
  primary_var_t* primary = polymec_malloc(sizeof(primary_var_t));
  primary->index = index;
  primary->size = size;
  primary->num_components = num_components;
  string_ptr_unordered_map_insert_with_kv_dtors(state->primary_map, (char*)var_name, primary, string_free, polymec_free);

  // Allocate or reallocate the solution vector.
  if (state->solution != NULL)
    polymec_free(state->solution);
  state->solution_size += size * num_components;
  state->solution = polymec_malloc(sizeof(real_t) * state->solution_size);
}

physics_state_t* physics_state_clone(physics_state_t* state)
{
  physics_state_t* state1 = physics_state_new();

  // Copy primary variable metadata.
  int pos = 0;
  char* var_name;
  void* entry;
  while (string_ptr_unordered_map_next(state->primary_map, &pos, &var_name, &entry))
  {
    primary_var_t* var = entry;
    physics_state_add_primary_with_index(state1, (const char*)var_name, var->index, 
                                         var->size, var->num_components);
  }

  // Copy primary variable data.
  ASSERT(state1->solution_size == state->solution_size);
  memcpy(state1->solution, state->solution, sizeof(real_t) * state1->solution_size);

  // Copy secondary variable metadata/data.
  pos = 0;
  while (string_ptr_unordered_map_next(state->secondary_map, &pos, &var_name, &entry))
  {
    secondary_var_t* var = entry;
    physics_state_add_secondary(state1, (const char*)var_name, 
                                var->size, var->num_components, var->kernel);
    secondary_var_t* var1 = *string_ptr_unordered_map_get(state1->secondary_map, var_name);
    memcpy(var1->data, var->data, sizeof(real_t) * var->size * var->num_components);
  }

  return state1;
}

real_t* physics_state_solution(physics_state_t* state)
{
  return state->solution;
}

void physics_state_add_primary(physics_state_t* state,
                               const char* var_name,
                               int size,
                               int num_components)
{
  physics_state_add_primary_with_index(state, var_name, state->max_index, size, num_components);
  ++(state->max_index);
}

void physics_state_extract_primary(physics_state_t* state, 
                                   const char* var_name, 
                                   int index,
                                   real_t* array)
{
  primary_var_t** var_p = (primary_var_t**)string_ptr_unordered_map_get(state->primary_map, (char*)var_name);
  if (var_p != NULL)
  {
    int var_index = (*var_p)->index;
    int num_comps = (*var_p)->num_components;
    for (int i = 0, j = 0; i < state->solution_size; i += var_index, ++j)
    {
      for (int c = 0; c < num_comps; ++c)
        array[index*j+c] = state->solution[i+c];
    }
  }
}

bool physics_state_has_primary(physics_state_t* state, const char* var_name)
{
  return (string_ptr_unordered_map_get(state->primary_map, (char*)var_name) != NULL);
}

bool physics_state_next_primary(physics_state_t* state,
                                int* pos,
                                char** var_name,
                                int* index,
                                int* size,
                                int* num_components)
{
  void* entry;
  bool result = string_ptr_unordered_map_next(state->primary_map, pos, var_name, &entry);
  if (result)
  {
    primary_var_t* var = entry;
    *index = var->index;
    *size = var->size;
    *num_components = var->num_components;
  }
  return result;
}

void physics_state_add_secondary(physics_state_t* state,
                                 const char* var_name,
                                 int size,
                                 int num_components,
                                 physics_kernel_t* kernel)
{
  ASSERT(size > 0);
  ASSERT(num_components > 0);
  ASSERT(kernel != NULL);
  secondary_var_t* secondary = polymec_malloc(sizeof(secondary_var_t));
  secondary->size = size;
  secondary->num_components = num_components;
  secondary->data = polymec_malloc(sizeof(real_t) * size * num_components);
  secondary->kernel = kernel;
  secondary->deps = string_array_new();
  string_ptr_unordered_map_insert_with_kv_dtors(state->secondary_map, (char*)var_name, secondary, string_free, secondary_var_free);
}

void physics_state_add_secondary_dep(physics_state_t* state,
                                     const char* var_name,
                                     const char* dep_name)
{
  secondary_var_t** var_p = (secondary_var_t**)string_ptr_unordered_map_get(state->secondary_map, (char*)var_name);
  ASSERT(var_p != NULL);
  secondary_var_t* var = *var_p;

#ifndef NDEBUG
  {
    // Make sure that this dependency doesn't appear already in the list of deps.
    char* deps[var->deps->size+1];
    for (int i = 0; i < var->deps->size; ++i)
      deps[i] = var->deps->data[i];
    deps[var->deps->size] = NULL;
    ASSERT(string_find_in_list(dep_name, (const char**)deps, false) == -1);
  }
#endif

  string_array_append(var->deps, (char*)dep_name);
}

void physics_state_extract_secondary(physics_state_t* state, 
                                     const char* var_name, 
                                     int index,
                                     real_t* array)
{
  secondary_var_t** var_p = (secondary_var_t**)string_ptr_unordered_map_get(state->secondary_map, (char*)var_name);
  if (var_p != NULL)
  {
    int size = (*var_p)->size;
    int num_comps = (*var_p)->num_components;
    real_t* data = (*var_p)->data;
    for (int i = 0; i < size; ++i)
    {
      for (int c = 0; c < num_comps; ++c)
        array[index*i+c] = data[num_comps*i+c];
    }
  }
}

bool physics_state_has_secondary(physics_state_t* state, const char* var_name)
{
  return (string_ptr_unordered_map_get(state->secondary_map, (char*)var_name) != NULL);
}

bool physics_state_next_secondary(physics_state_t* state,
                                  int* pos,
                                  char** var_name,
                                  int* size,
                                  int* num_components,
                                  physics_kernel_t** kernel)
{
  void* entry;
  bool result = string_ptr_unordered_map_next(state->secondary_map, pos, var_name, &entry);
  if (result)
  {
    secondary_var_t* var = entry;
    *size = var->size;
    *num_components = var->num_components;
    *kernel = var->kernel;
  }
  return result;
}

bool physics_state_next_secondary_dep(physics_state_t* state,
                                      char* var_name,
                                      int* pos,
                                      char** dep_name)
{
  secondary_var_t** var_p = (secondary_var_t**)string_ptr_unordered_map_get(state->secondary_map, (char*)var_name);
  if (var_p == NULL) 
    return false;

  secondary_var_t* var = *var_p;
  if (*pos < var->deps->size)
  {
    *dep_name = var->deps->data[*pos];
    ++(*pos);
    return true;
  }
  else
    return false;
}

