// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/array.h"
#include "model/model_probe.h"

typedef void (*callback_t)(void*, real_t, int, size_t*, real_t*);
DEFINE_ARRAY(callback_array, callback_t)

struct model_probe_t
{
  char* name;
  int data_rank;
  size_t* data_shape;
  void* context;
  model_probe_vtable vtable;

  ptr_array_t* callback_contexts;
  callback_array_t* callbacks;
};

model_probe_t* model_probe_new(const char* name, 
                               int data_rank,
                               size_t* data_shape,
                               void* context, 
                               model_probe_vtable vtable)
{
  ASSERT(data_rank >= 0);
  ASSERT((data_shape != NULL) || (data_rank == 0));
  ASSERT(vtable.acquire != NULL);

  model_probe_t* probe = polymec_malloc(sizeof(model_probe_t));
  probe->name = string_dup(name);
  probe->data_rank = data_rank;
  probe->data_shape = polymec_malloc(sizeof(size_t) * data_rank);
  probe->context = context;
  for (int i = 0; i < data_rank; ++i)
  {
    ASSERT(probe->data_shape[i] > 0);
    probe->data_shape[i] = data_shape[i];
  }
  probe->vtable = vtable;
  probe->callback_contexts = ptr_array_new();
  probe->callbacks = callback_array_new();
  return probe;
}

void model_probe_free(model_probe_t* probe)
{
  string_free(probe->name);
  if ((probe->context != NULL) && (probe->vtable.dtor != NULL))
    probe->vtable.dtor(probe->context);
  polymec_free(probe->data_shape);
  ptr_array_free(probe->callback_contexts);
  callback_array_free(probe->callbacks);
  polymec_free(probe);
}

char* model_probe_name(model_probe_t* probe)
{
  return probe->name;
}

int model_probe_data_rank(model_probe_t* probe)
{
  return probe->data_rank;
}

size_t* model_probe_data_shape(model_probe_t* probe)
{
  return probe->data_shape;
}

real_t* model_probe_new_array(model_probe_t* probe)
{
  size_t n = 1;
  for (int i = 0; i < probe->data_rank; ++i)
    n *= probe->data_shape[i];
  return polymec_malloc(sizeof(real_t) * n);
}

void model_probe_acquire(model_probe_t* probe, real_t t, real_t* data)
{
  probe->vtable.acquire(probe->context, t, probe->data_rank, probe->data_shape, data);
  for (size_t i = 0; i < probe->callbacks->size; ++i)
    probe->callbacks->data[i](probe->callback_contexts->data[i], t, probe->data_rank, probe->data_shape, data);
}

void model_probe_add_callback(model_probe_t* probe, 
                              void* context, 
                              void (*callback)(void* context,
                                               real_t t, 
                                               int rank,
                                               size_t* shape,
                                               real_t* data))
{
  ptr_array_append(probe->callback_contexts, context);
  callback_array_append(probe->callbacks, callback);
}
