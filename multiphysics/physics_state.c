// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/timer.h"
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

  // Sorted list of secondary variables to update.
  string_array_t* sorted_secondaries;
  bool updating_dependencies;
};

physics_state_t* physics_state_new()
{
  physics_state_t* state = polymec_malloc(sizeof(physics_state_t));
  state->solution = NULL;
  state->solution_size = 0;
  state->max_index = 0;
  state->primary_map = string_ptr_unordered_map_new();
  state->secondary_map = string_ptr_unordered_map_new();
  state->sorted_secondaries = NULL;
  state->updating_dependencies = false;
  return state;
}

void physics_state_free(physics_state_t* state)
{
  if (state->solution != NULL)
    polymec_free(state->solution);
  string_ptr_unordered_map_free(state->primary_map);
  string_ptr_unordered_map_free(state->secondary_map);
  if (state->sorted_secondaries != NULL)
    string_array_free(state->sorted_secondaries);
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
  // This function must not be called after we have extracted data from 
  // the state.
  ASSERT(state->sorted_secondaries == NULL); 

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
  // This function must not be called after we have extracted data from 
  // the state.
  ASSERT(state->sorted_secondaries == NULL); 

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
  // This function must not be called after we have extracted data from 
  // the state.
  ASSERT(state->sorted_secondaries == NULL); 

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

static void update_dependencies(physics_state_t* state, real_t t)
{
  START_FUNCTION_TIMER();
  // Generate a sorted list of our secondary variables if needed.
  if (state->sorted_secondaries == NULL)
  {
    // Make a list of all the secondary variables, and index them.
    string_array_t* unsorted_vars = string_array_new();
    string_int_unordered_map_t* var_indices = string_int_unordered_map_new();
    int pos = 0;
    char* var;
    void* entry;
    while (string_ptr_unordered_map_next(state->secondary_map, &pos, &var, &entry))
    {
      string_int_unordered_map_insert(var_indices, var, unsorted_vars->size);
      string_array_append(unsorted_vars, var);
    }

    // Create an adjacency graph representing the dependencies. Note that in 
    // canonical graph theory, vertex "v depends on w" means w -> v.
    int num_secondaries = unsorted_vars->size;
    adj_graph_t* graph = adj_graph_new(MPI_COMM_SELF, num_secondaries);
    string_array_t* deps = string_array_new();
    for (int v = 0; v < unsorted_vars->size; ++v)
    {
      string_array_clear(deps);
      int pos = 0; 
      char *dep;
      while (physics_state_next_secondary_dep(state, unsorted_vars->data[v], &pos, &dep))
        string_array_append(deps, dep);
      for (int j = 0; j < deps->size; ++j)
      {
        int w = *string_int_unordered_map_get(var_indices, deps->data[j]);
        int k = adj_graph_num_edges(graph, w);
        adj_graph_set_num_edges(graph, w, k + 1);
        int* edges = adj_graph_edges(graph, w);
        edges[k] = v;
      }
    }
    string_array_free(deps);
    string_int_unordered_map_free(var_indices);

    // Sort the graph.
    int sorted_indices[num_secondaries];
    bool sorted = adj_graph_sort(graph, sorted_indices);
    if (!sorted)
      polymec_error("physics_state: cycle detected in dependencies for secondary variables.");

    // Assemble a proper list of secondaries in sorted order.
    state->sorted_secondaries = string_array_new();
    for (int i = 0; i < num_secondaries; ++i)
      string_array_append(state->sorted_secondaries, unsorted_vars->data[sorted_indices[i]]);

    string_array_free(unsorted_vars);
  }

  // Now update each variable.
  state->updating_dependencies = true;
  int pos = 0, size, num_components;
  char* var_name;
  physics_kernel_t* kernel;
  while (physics_state_next_secondary(state, &pos, &var_name, &size, &num_components, &kernel))
    physics_kernel_update_secondary(kernel, var_name, t, state);
  state->updating_dependencies = false;
  STOP_FUNCTION_TIMER();
}

void physics_state_extract_secondary(physics_state_t* state, 
                                     const char* var_name, 
                                     real_t t,
                                     int index,
                                     real_t* array)
{
  // Make sure the dependencies are sorted and updated, but only if we're 
  // not already being invoked during a call to update_dependencies.
  if (!state->updating_dependencies)
    update_dependencies(state, t);

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
                                      const char* var_name,
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

real_t* physics_state_secondary(physics_state_t* state,
                                const char* var_name)
{
  secondary_var_t** var_p = (secondary_var_t**)string_ptr_unordered_map_get(state->secondary_map, (char*)var_name);
  if (var_p == NULL) 
    return NULL;

  return (*var_p)->data;
}

