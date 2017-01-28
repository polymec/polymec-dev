// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "model/model_probe.h"

struct model_datum_t
{
  int rank;
  size_t* shape;
  real_t* data;
};

model_datum_t* model_datum_new(int rank, size_t* shape)
{
  ASSERT(rank >= 0);
  ASSERT((shape != NULL) || (rank == 0));

  model_datum_t* datum = polymec_malloc(sizeof(model_datum_t));
  datum->rank = rank;
  datum->shape = polymec_malloc(sizeof(size_t) * rank);
  size_t size = 1;
  for (int i = 0; i < rank; ++i)
  {
    ASSERT(datum->shape[i] > 0);
    datum->shape[i] = shape[i];
    size *= shape[i];
  }
  datum->data = polymec_malloc(sizeof(real_t) * size);
  return datum;
}

void model_datum_free(model_datum_t* datum)
{
  polymec_free(datum->data);
  polymec_free(datum->shape);
  polymec_free(datum);
}

int model_datum_rank(model_datum_t* datum)
{
  return datum->rank;
}

size_t* model_datum_shape(model_datum_t* datum)
{
  return datum->shape;
}

real_t* model_datum_data(model_datum_t* datum)
{
  return datum->data;
}

struct model_probe_t
{
  char* name;
  int datum_rank;
  size_t* datum_shape;
  void* context;
  model_probe_vtable vtable;
};

model_probe_t* model_probe_new(const char* name, 
                               int datum_rank,
                               size_t* datum_shape,
                               void* context, 
                               model_probe_vtable vtable)
{
  ASSERT(datum_rank >= 0);
  ASSERT((datum_shape != NULL) || (datum_rank == 0));
  ASSERT(vtable.acquire != NULL);

  model_probe_t* probe = polymec_malloc(sizeof(model_probe_t));
  probe->name = string_dup(name);
  probe->datum_rank = datum_rank;
  probe->datum_shape = polymec_malloc(sizeof(size_t) * datum_rank);
  for (int i = 0; i < datum_rank; ++i)
  {
    ASSERT(probe->datum_shape[i] > 0);
    probe->datum_shape[i] = datum_shape[i];
  }
  probe->vtable = vtable;
  return probe;
}

void model_probe_free(model_probe_t* probe)
{
  string_free(probe->name);
  if ((probe->context != NULL) && (probe->vtable.dtor != NULL))
    probe->vtable.dtor(probe->context);
  polymec_free(probe->datum_shape);
  polymec_free(probe);
}

char* model_probe_name(model_probe_t* probe)
{
  return probe->name;
}

model_datum_t* model_probe_new_datum(model_probe_t* probe)
{
  return model_datum_new(probe->datum_rank, probe->datum_shape);
}

void model_probe_acquire(model_probe_t* probe, real_t t, model_datum_t* datum)
{
  probe->vtable.acquire(probe->context, t, datum);
}

