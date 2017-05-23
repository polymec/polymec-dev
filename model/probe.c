// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/array.h"
#include "model/probe.h"

probe_data_t* probe_data_new(int rank, size_t* shape)
{
  ASSERT(rank >= 0);
  probe_data_t* data = polymec_malloc(sizeof(probe_data_t));
  data->rank = rank;
  data->shape = polymec_malloc(sizeof(size_t) * rank);
  size_t n = 1;
  for (int i = 0; i < rank; ++i)
  {
    data->shape[i] = shape[i];
    n *= shape[i];
  }
  data->data = polymec_malloc(sizeof(real_t) * n);
  return data;
}

void probe_data_free(probe_data_t* data)
{
  polymec_free(data->data);
  polymec_free(data->shape);
  polymec_free(data);
}

struct probe_t
{
  char* name;
  char* data_name;
  int rank;
  size_t* shape;
  void* context;
  probe_vtable vtable;
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
  return probe;
}

void probe_free(probe_t* probe)
{
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

probe_data_t* probe_acquire(probe_t* probe, real_t t)
{
  probe_data_t* data = probe_data_new(probe->rank, probe->shape);
  probe->vtable.acquire(probe->context, t, data);
  data->time = t;
  return data;
}

