// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/array.h"
#include "model/probe.h"
#include "model/model.h"

probe_data_t* probe_data_new(int rank, size_t* shape)
{
  ASSERT(rank >= 0);
  probe_data_t* data = polymec_malloc(sizeof(probe_data_t));
  data->rank = rank;
  data->shape = polymec_malloc(sizeof(size_t) * rank);
  for (int i = 0; i < data->rank; ++i)
    data->shape[i] = shape[i];
  size_t data_size = probe_data_size(data);
  data->data = polymec_malloc(sizeof(real_t) * data_size);
  return data;
}

size_t probe_data_size(probe_data_t* data)
{
  size_t n = 1;
  for (int i = 0; i < data->rank; ++i)
    n *= data->shape[i];
  return n;
}

void probe_data_free(probe_data_t* data)
{
  polymec_free(data->data);
  polymec_free(data->shape);
  polymec_free(data);
}

// Here's an acquisition callback.
typedef struct
{
  void* context;
  void (*function)(void* context, real_t t, probe_data_t* data);
  void (*dtor)(void* context);
} acq_callback_t;

static void free_callback_context(acq_callback_t callback)
{
  if ((callback.context != NULL) && (callback.dtor != NULL))
    callback.dtor(callback.context);
}

DEFINE_ARRAY(acq_callback_array, acq_callback_t)

struct probe_t
{
  char* name;
  char* data_name;
  int rank;
  size_t* shape;
  void* context;
  probe_vtable vtable;
  acq_callback_array_t* callbacks;
};

probe_t* probe_new(const char* name,
                   const char* data_name,
                   int rank,
                   size_t* shape,
                   void* context,
                   probe_vtable vtable)
{
  ASSERT(rank >= 0);
  ASSERT((shape != NULL) || (rank == 0));
  ASSERT(vtable.acquire != NULL);

  probe_t* probe = polymec_malloc(sizeof(probe_t));
  probe->name = string_dup(name);
  probe->data_name = string_dup(data_name);
  probe->rank = rank;
  probe->shape = polymec_malloc(sizeof(size_t) * rank);
  probe->context = context;
  for (int i = 0; i < rank; ++i)
  {
    ASSERT(shape[i] > 0);
    probe->shape[i] = shape[i];
  }
  probe->vtable = vtable;
  probe->callbacks = acq_callback_array_new();
  return probe;
}

void probe_free(probe_t* probe)
{
  acq_callback_array_free(probe->callbacks);
  string_free(probe->name);
  string_free(probe->data_name);
  if ((probe->context != NULL) && (probe->vtable.dtor != NULL))
    probe->vtable.dtor(probe->context);
  polymec_free(probe->shape);
  polymec_free(probe);
}

char* probe_name(probe_t* probe)
{
  return probe->name;
}

char* probe_data_name(probe_t* probe)
{
  return probe->data_name;
}

void* probe_context(probe_t* probe)
{
  return probe->context;
}

probe_data_t* probe_acquire(probe_t* probe, real_t t)
{
  // Allocate storage for and acquire the data.
  probe_data_t* data = probe_data_new(probe->rank, probe->shape);
  probe->vtable.acquire(probe->context, t, data);
  data->time = t;

  // Call any callbacks we have with the time and the data.
  for (size_t i = 0; i < probe->callbacks->size; ++i)
  {
    acq_callback_t* callback = &(probe->callbacks->data[i]);
    callback->function(callback->context, t, data);
  }

  // Return the data.
  return data;
}

void probe_set_model(probe_t* probe, model_t* model);
void probe_set_model(probe_t* probe, model_t* model)
{
  if (probe->vtable.set_model != NULL)
    probe->vtable.set_model(probe_context(probe), model_context(model));
}

void probe_postprocess(probe_t* probe, real_array_t* times, probe_data_array_t* data)
{
  if ((probe->vtable.postprocess != NULL) && (data != NULL))
    probe->vtable.postprocess(probe->context, times, data);
}

void probe_on_acquire(probe_t* probe,
                      void* context,
                      void (*function)(void* context, real_t t, probe_data_t* data),
                      void (*dtor)(void* context))
{
  acq_callback_t callback = {.context = context, .function = function, .dtor = dtor};
  acq_callback_array_append_with_dtor(probe->callbacks, callback, free_callback_context);
}

