// Copyright 2012-2013 Jeffrey Johnson.
// 
// This file is part of Polymec, and is licensed under the Apache License, 
// Version 2.0 (the "License"); you may not use this file except in 
// compliance with the License. You may may find the text of the license in 
// the LICENSE file at the top-level source directory, or obtain a copy of 
// it at
// 
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

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

