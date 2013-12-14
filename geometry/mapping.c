// Copyright (c) 2012-2013, Jeffrey N. Johnson
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this 
// list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice, 
// this list of conditions and the following disclaimer in the documentation 
// and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include <gc/gc.h>
#include "geometry/mapping.h"

struct mapping_t 
{
  char* name;
  void* context;
  mapping_vtable vtable;
};

static void mapping_free(void* ctx, void* dummy)
{
  mapping_t* mapping = (mapping_t*)ctx;
  if (mapping->vtable.dtor)
    free(mapping->context);
  free(mapping->name);
}

mapping_t* mapping_new(const char* name, void* context, mapping_vtable vtable)
{
  ASSERT(context != NULL);
  ASSERT(vtable.map != NULL);
  ASSERT(vtable.jacobian != NULL);
  mapping_t* m = GC_MALLOC(sizeof(mapping_t));
  m->name = string_dup(name);
  m->context = context;
  m->vtable = vtable;
  GC_register_finalizer(m, mapping_free, m, NULL, NULL);
  return m;
}

const char* mapping_name(mapping_t* mapping)
{
  return mapping->name;
}

void* mapping_context(mapping_t* mapping)
{
  return mapping->context;
}

void mapping_map(mapping_t* mapping, point_t* x, point_t* y)
{
  mapping->vtable.map(mapping->context, x, y);
}

void mapping_compute_jacobian(mapping_t* mapping, point_t* x, double* J)
{
  mapping->vtable.jacobian(mapping->context, x, J);
}

