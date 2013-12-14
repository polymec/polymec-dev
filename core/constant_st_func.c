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

