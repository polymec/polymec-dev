// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "gc/gc.h"
#include "model/aux_state.h"

struct aux_state_t 
{
  int_array_t* types;
  int_array_t* indices;
  int_array_t* offsets;
  real_array_t* values;
};

static void aux_state_free(void* ctx, void* dummy)
{
  aux_state_t* state = ctx;
  int_array_free(state->types);
  int_array_free(state->indices);
  int_array_free(state->offsets);
  real_array_free(state->values);
}

aux_state_t* aux_state_new()
{
  aux_state_t* state = GC_MALLOC(sizeof(aux_state_t));
  state->types = int_array_new();
  state->indices = int_array_new();
  state->offsets = int_array_new();
  int_array_append(state->offsets, 0);
  state->values = real_array_new();
  GC_register_finalizer(state, aux_state_free, state, NULL, NULL);
  return state;
}

aux_state_t* aux_state_clone(aux_state_t* state)
{
  aux_state_t* clone = aux_state_new();
  int_array_resize(clone->types, state->types->size);
  memcpy(clone->types->data, state->types->data, sizeof(int) * state->types->size);
  int_array_resize(clone->indices, state->indices->size);
  memcpy(clone->indices->data, state->indices->data, sizeof(int) * state->indices->size);
  int_array_resize(clone->offsets, state->offsets->size);
  memcpy(clone->offsets->data, state->offsets->data, sizeof(int) * state->offsets->size);
  real_array_resize(clone->values, state->values->size);
  memcpy(clone->values->data, state->values->data, sizeof(real_t) * state->values->size);
  return clone;
}

int aux_state_size(aux_state_t* state)
{
  return state->types->size;
}

// Returns the variable offset for the given index, or -1 if it is not found.
static inline int index_map(aux_state_t* state, int index)
{
  for (int i = 0; i < state->indices->size; ++i)
  {
    if (state->indices->data[i] == index)
      return i;
  }
  return -1;
}

static inline bool aux_state_index_found(aux_state_t* state, int index)
{
  return (index_map(state, index) != -1);
}

void aux_state_add_scalar(aux_state_t* state, int index)
{
  ASSERT(!aux_state_index_found(state, index));
  int_array_append(state->types, AUX_STATE_SCALAR);
  int_array_append(state->indices, index);
  int_array_append(state->offsets, state->offsets->data[state->offsets->size-1] + 1);
  real_array_append(state->values, 0.0);
}

void aux_state_add_vector(aux_state_t* state, int index)
{
  ASSERT(!aux_state_index_found(state, index));
  int_array_append(state->types, AUX_STATE_VECTOR);
  int_array_append(state->indices, index);
  int_array_append(state->offsets, state->offsets->data[state->offsets->size-1] + 3);
  real_array_append(state->values, 0.0);
  real_array_append(state->values, 0.0);
  real_array_append(state->values, 0.0);
}

void aux_state_add_tensor2(aux_state_t* state, int index)
{
  ASSERT(!aux_state_index_found(state, index));
  int_array_append(state->types, AUX_STATE_TENSOR2);
  int_array_append(state->indices, index);
  int_array_append(state->offsets, state->offsets->data[state->offsets->size-1] + 9);
  real_array_append(state->values, 0.0);
  real_array_append(state->values, 0.0);
  real_array_append(state->values, 0.0);
  real_array_append(state->values, 0.0);
  real_array_append(state->values, 0.0);
  real_array_append(state->values, 0.0);
  real_array_append(state->values, 0.0);
  real_array_append(state->values, 0.0);
  real_array_append(state->values, 0.0);
}

void aux_state_add_sym_tensor2(aux_state_t* state, int index)
{
  ASSERT(!aux_state_index_found(state, index));
  int_array_append(state->types, AUX_STATE_SYM_TENSOR2);
  int_array_append(state->indices, index);
  int_array_append(state->offsets, state->offsets->data[state->offsets->size-1] + 6);
  real_array_append(state->values, 0.0);
  real_array_append(state->values, 0.0);
  real_array_append(state->values, 0.0);
  real_array_append(state->values, 0.0);
  real_array_append(state->values, 0.0);
  real_array_append(state->values, 0.0);
}

void aux_state_set_scalar(aux_state_t* state, int index, real_t value)
{
  int I = index_map(state, index);
  ASSERT(I != -1);
  int offset = state->offsets->data[I];
  state->values->data[offset] = value;
}

real_t aux_state_scalar(aux_state_t* state, int index)
{
  int I = index_map(state, index);
  ASSERT(I != -1);
  ASSERT(state->types->data[I] == AUX_STATE_SCALAR);
  int offset = state->offsets->data[I];
  return state->values->data[offset];
}

void aux_state_set_vector(aux_state_t* state, int index, vector_t* value)
{
  int I = index_map(state, index);
  ASSERT(I != -1);
  int offset = state->offsets->data[I];
  state->values->data[offset]   = value->x;
  state->values->data[offset+1] = value->y;
  state->values->data[offset+2] = value->z;
}

vector_t* aux_state_vector(aux_state_t* state, int index)
{
  int I = index_map(state, index);
  ASSERT(I != -1);
  ASSERT(state->types->data[I] == AUX_STATE_VECTOR);
  int offset = state->offsets->data[I];
  return (vector_t*)&(state->values->data[offset]);
}

void aux_state_set_tensor2(aux_state_t* state, int index, tensor2_t* value)
{
  int I = index_map(state, index);
  ASSERT(I != -1);
  int offset = state->offsets->data[I];
  memcpy(&(state->values->data[offset]), (real_t*)value, 9 * sizeof(real_t));
}

tensor2_t* aux_state_tensor2(aux_state_t* state, int index)
{
  int I = index_map(state, index);
  ASSERT(I != -1);
  ASSERT(state->types->data[I] == AUX_STATE_TENSOR2);
  int offset = state->offsets->data[I];
  return (tensor2_t*)&(state->values->data[offset]);
}

void aux_state_set_sym_tensor2(aux_state_t* state, int index, sym_tensor2_t* value)
{
  int I = index_map(state, index);
  ASSERT(I != -1);
  int offset = state->offsets->data[I];
  memcpy(&(state->values->data[offset]), (real_t*)value, 6 * sizeof(real_t));
}

sym_tensor2_t* aux_state_sym_tensor2(aux_state_t* state, int index)
{
  int I = index_map(state, index);
  ASSERT(I != -1);
  ASSERT(state->types->data[I] == AUX_STATE_SYM_TENSOR2);
  int offset = state->offsets->data[I];
  return (sym_tensor2_t*)&(state->values->data[offset]);
}

aux_state_var_t aux_state_type(aux_state_t* state, int index)
{
  int I = index_map(state, index);
  ASSERT(I != -1);
  return (aux_state_var_t)state->types->data[I];
}

bool aux_state_has_scalar(aux_state_t* state, int index)
{
  int I = index_map(state, index);
  if (I == -1) 
    return false;
  else
    return (state->types->data[I] == AUX_STATE_SCALAR);
}

bool aux_state_has_vector(aux_state_t* state, int index)
{
  int I = index_map(state, index);
  if (I == -1) 
    return false;
  else
    return (state->types->data[I] == AUX_STATE_VECTOR);
}

bool aux_state_has_tensor2(aux_state_t* state, int index)
{
  int I = index_map(state, index);
  if (I == -1) 
    return false;
  else
    return (state->types->data[I] == AUX_STATE_TENSOR2);
}

bool aux_state_has_sym_tensor2(aux_state_t* state, int index)
{
  int I = index_map(state, index);
  if (I == -1) 
    return false;
  else
    return (state->types->data[I] == AUX_STATE_SYM_TENSOR2);
}

