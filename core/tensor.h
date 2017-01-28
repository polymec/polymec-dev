// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_TENSOR_H
#define POLYMEC_TENSOR_H

#include "core/polymec.h"

// A real-valued general tensor of specifiable rank in a space of specifiable 
// dimension. The "shape" of a tensor is an array whose components are the 
// dimensions of its indices. This is a tensor in the broad mathematical sense, 
// not a spatial one. These tend to be larger than the rank-specific spatial 
// tensors like tensor2s.
typedef struct tensor_t tensor_t;

// Allocates a tensor with the given rank and shape.
tensor_t* tensor_new(int rank, size_t* shape);

// Returns a deep copy of the given tensor.
tensor_t* tensor_clone(tensor_t* t);

// Destroys the given tensor.
void tensor_free(tensor_t* t);

// Returns the rank of the given tensor.
int tensor_rank(tensor_t* t);

// Returns an internal pointer to an array describing the shape of the 
// given tensor.
size_t* tensor_shape(tensor_t* shape);

// Returns an internal pointer to the data within the given tensor.
real_t* tensor_data(tensor_t* t);

// Returns the total number of components in the tensor.
size_t tensor_size(tensor_t* t);

#endif

