// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_AUX_STATE_H
#define POLYMEC_AUX_STATE_H

#include "core/point.h"
#include "core/tensor2.h"

// An auxiliary state is an object that stores information pertaining to a 
// state at a fixed spatial location. The nature of the state depends on the 
// particular algorithm using it. Really, it's just an opaque array of 
// numerical values that are indexed and retrieved in a manner defined by the 
// model using it. Objects of this type are garbage-collected.
typedef struct aux_state_t aux_state_t;

// The following enumerated type provides types of quantities that can be
// stored in an auxiliary state.
typedef enum
{
  AUX_STATE_SCALAR,
  AUX_STATE_VECTOR,
  AUX_STATE_TENSOR2,
  AUX_STATE_SYM_TENSOR2
} aux_state_var_t;

// Creates a new empty auxiliary state.
aux_state_t* aux_state_new(void);

// Creates a (deep) copy of the given auxiliary state.
aux_state_t* aux_state_clone(aux_state_t* state);

// Returns the number of quantities stored in this auxiliary state.
int aux_state_size(aux_state_t* state);

// Appends storage for a scalar-valued quantity to the state, returning the 
// index at which it is stored.
void aux_state_add_scalar(aux_state_t* state, int index);

// Appends storage for a vector-valued quantity to the state, returning the 
// index at which it is stored.
void aux_state_add_vector(aux_state_t* state, int index);

// Appends storage for a (rank 2) tensor-valued quantity to the state, 
// returning the index at which it is stored.
void aux_state_add_tensor2(aux_state_t* state, int index);

// Appends storage for a (rank 2) symmetric-tensor quantity to the state, 
// returning the index at which it is stored.
void aux_state_add_sym_tensor2(aux_state_t* state, int index);

// Sets the value of the scalar at the given index within the state.
void aux_state_set_scalar(aux_state_t* state, int index, real_t value);

// Returns the scalar value at the given index within the state.
real_t aux_state_scalar(aux_state_t* state, int index);

// Sets the value of the vector at the given index within the state.
void aux_state_set_vector(aux_state_t* state, int index, vector_t* value);

// Returns an internal pointer to the vector value at the given index within 
// the state.
vector_t* aux_state_vector(aux_state_t* state, int index);

// Sets the value of the (rank-2) tensor at the given index within the state.
void aux_state_set_tensor2(aux_state_t* state, int index, tensor2_t* value);

// Returns an internal pointer to the rank 2 tensor value at the given index 
// within the state.
tensor2_t* aux_state_tensor2(aux_state_t* state, int index);

// Sets the value of the (rank-2) symmetric tensor at the given index within 
// the state.
void aux_state_set_sym_tensor2(aux_state_t* state, int index, sym_tensor2_t* value);

// Returns an internal pointer to the symmetric 2 tensor value at the given 
// index within the state.
sym_tensor2_t* aux_state_sym_tensor2(aux_state_t* state, int index);

// Returns the type of quantity stored in the given index within the state.
aux_state_var_t aux_state_type(aux_state_t* state, int index);

// Returns true if the quantity at the given index is a scalar, false if not.
bool aux_state_has_scalar(aux_state_t* state, int index);

// Returns true if the quantity at the given index is a vector, false if not.
bool aux_state_has_vector(aux_state_t* state, int index);

// Returns true if the quantity at the given index is a rank 2 tensor, false if not.
bool aux_state_has_tensor2(aux_state_t* state, int index);

// Returns true if the quantity at the given index is a rank 2 symmetric tensor, 
// false if not.
bool aux_state_has_sym_tensor2(aux_state_t* state, int index);

#endif
