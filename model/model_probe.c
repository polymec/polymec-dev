// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "model/model_probe.h"

struct model_probe_t
{
  char* name;
  size_t size;
  void* context;
  model_probe_vtable vtable;
};

model_probe_t* model_probe_new(const char* name, 
                               size_t datum_size,
                               void* context, 
                               model_probe_vtable vtable)
{
  ASSERT(datum_size > 0);
  ASSERT(vtable.acquire != NULL);

  model_probe_t* probe = polymec_malloc(sizeof(model_probe_t));
  probe->name = string_dup(name);
  probe->size = datum_size;
  probe->vtable = vtable;
  return probe;
}

void model_probe_free(model_probe_t* probe)
{
  string_free(probe->name);
  if ((probe->context != NULL) && (probe->vtable.dtor != NULL))
    probe->vtable.dtor(probe->context);
  polymec_free(probe);
}

char* model_probe_name(model_probe_t* probe)
{
  return probe->name;
}

size_t model_probe_datum_size(model_probe_t* probe)
{
  return probe->size;
}

void model_probe_acquire(model_probe_t* probe, real_t t, real_t* datum)
{
  probe->vtable.acquire(probe->context, t, datum);
}

