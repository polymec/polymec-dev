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

#include "integrators/integrator.h"

struct integrator_t 
{
  char* name;
  void* context;
  integrator_vtable vtable;
  int order;
  integrator_type_t type;
  int N; // Number of unknowns.
};

integrator_t* integrator_new(const char* name, 
                             void* context,
                             integrator_vtable vtable,
                             int order,
                             integrator_type_t type,
                             int N)
{
  ASSERT(vtable.step != NULL);
  ASSERT(order > 0);
  ASSERT(N > 0);
  integrator_t* integ = malloc(sizeof(integrator_t));
  integ->name = string_dup(name);
  integ->context = context;
  integ->vtable = vtable;
  integ->order = order;
  integ->type = type;
  integ->N = N;
  return integ;
}

void integrator_free(integrator_t* integrator)
{
  if ((integrator->context != NULL) && (integrator->vtable.dtor != NULL))
    integrator->vtable.dtor(integrator->context);
  free(integrator->name);
  free(integrator);
}

char* integrator_name(integrator_t* integrator)
{
  return integrator->name;
}

void* integrator_context(integrator_t* integrator)
{
  return integrator->context;
}

int integrator_order(integrator_t* integrator)
{
  return integrator->order;
}

integrator_type_t integrator_type(integrator_t* integrator)
{
  return integrator->type;
}

int integrator_N(integrator_t* integrator)
{
  return integrator->N;
}

void integrator_init(integrator_t* integrator, double t, double* X0)
{
  integrator->vtable.init(integrator->context, t, integrator->N, X0);
}

void integrator_step(integrator_t* integrator, double t1, double t2, double* X)
{
  ASSERT(t2 > t1);
  integrator->vtable.step(integrator->context, t1, t2, integrator->N, X);
}

