// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/tensor.h"

struct tensor_t
{
  int rank;
  size_t* shape;
  real_t* data;
};

tensor_t* tensor_new(int rank, size_t* shape)
{
  ASSERT(rank >= 0);
  ASSERT((shape != NULL) || (rank == 0));

  tensor_t* t = polymec_malloc(sizeof(tensor_t));
  t->rank = rank;
  t->shape = polymec_malloc(sizeof(size_t) * rank);
  size_t size = 1;
  for (int i = 0; i < rank; ++i)
  {
    ASSERT(shape[i] > 0);
    t->shape[i] = shape[i];
    size *= shape[i];
  }
  t->data = polymec_malloc(sizeof(real_t) * size);
  return t;
}

tensor_t* tensor_clone(tensor_t* t)
{
  tensor_t* clone = tensor_new(t->rank, t->shape);
  size_t s = tensor_size(t);
  memcpy(clone->data, t->data, sizeof(real_t) * s);
  return clone;
}

void tensor_free(tensor_t* t)
{
  polymec_free(t->data);
  polymec_free(t->shape);
  polymec_free(t);
}

int tensor_rank(tensor_t* t)
{
  return t->rank;
}

size_t* tensor_shape(tensor_t* t)
{
  return t->shape;
}

real_t* tensor_data(tensor_t* t)
{
  return t->data;
}

size_t tensor_size(tensor_t* t)
{
  size_t s = 1;
  for (int i = 0; i < t->rank; ++i)
    s *= t->shape[i];
  return s;
}
