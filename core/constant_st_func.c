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

#include <stdlib.h>
#include "core/constant_st_func.h"

typedef struct
{
  int num_comp;
  double *comp;
} const_st_func_t;

static void constant_eval(void* ctx, point_t* x, double t, double* res)
{
  const_st_func_t* f = (const_st_func_t*)ctx;
  for (int i = 0; i < f->num_comp; ++i)
    res[i] = f->comp[i];
}

static void constant_dtor(void* ctx)
{
  const_st_func_t* f = (const_st_func_t*)ctx;
  free(f->comp);
  free(f);
}

static st_func_t* create_constant_st_func(int num_comp, double comp[])
{
  st_vtable vtable = {.eval = constant_eval, .dtor = constant_dtor};
  char name[1024];
  snprintf(name, 1024, "constant space-time function"); // FIXME
  const_st_func_t* f = malloc(sizeof(const_st_func_t));
  f->num_comp = num_comp;
  f->comp = malloc(num_comp*sizeof(double));
  for (int i = 0; i < num_comp; ++i)
    f->comp[i] = comp[i];
  return st_func_new(name, (void*)f, vtable, ST_HOMOGENEOUS, ST_CONSTANT, num_comp);
}

st_func_t* constant_st_func_new(int num_comp, double comp[])
{
  st_func_t* func = create_constant_st_func(num_comp, comp);

  // Just to be complete, we register the 3-component zero function as this 
  // function's gradient.
  double zeros[3*num_comp];
  memset(zeros, 0, 3*num_comp*sizeof(double));
  st_func_t* zero = create_constant_st_func(3*num_comp, zeros);
  st_func_register_deriv(func, 1, zero);

  return func;
}

// For free!
sp_func_t* constant_sp_func_new(int num_comp, double comp[])
{
  return st_func_freeze(constant_st_func_new(num_comp, comp), 0.0);
}

