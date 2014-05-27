// Copyright (c) 2012-2014, Jeffrey N. Johnson
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

#include "geometry/difference.h"

typedef struct
{
  sp_func_t *s1, *s2;
} diff_t;

static void diff_eval(void* ctx, point_t* x, real_t* result)
{
  diff_t* diff = ctx;
  real_t v1, v2;
  sp_func_eval(diff->s1, x, &v1);
  sp_func_eval(diff->s2, x, &v2);
  result[0] = MAX(v1, -v2);
}

static void diff_eval_gradient(void* ctx, point_t* x, real_t* result)
{
  diff_t* diff = ctx;
  real_t v1[3], v2[3];
  sp_func_eval(diff->s1, x, v1);
  real_t v1mag2 = v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2];
  sp_func_eval(diff->s2, x, v2);
  real_t v2mag2 = v2[0]*v2[0] + v2[1]*v2[1] + v2[2]*v2[2];

  if (v1mag2 > v2mag2)
  {
    result[0] = v1[0];
    result[1] = v1[1];
    result[2] = v1[2];
  }
  else
  {
    result[0] = v2[0];
    result[1] = v2[1];
    result[2] = v2[2];
  }
}

sp_func_t* difference_new(sp_func_t* surface1, sp_func_t* surface2)
{
  ASSERT(surface1 != NULL);
  ASSERT(sp_func_num_comp(surface1) == 1);
  ASSERT(surface2 != NULL);
  ASSERT(sp_func_num_comp(surface2) == 1);

  diff_t* diff = polymec_malloc(sizeof(diff_t));
  diff->s1 = surface1;
  diff->s2 = surface2;

  char diff_str[1024];
  sprintf(diff_str, "Difference"); // FIXME: Not very helpful.
  sp_vtable vtable = {.eval = diff_eval, .dtor = polymec_free};
  sp_func_t* difference = sp_func_new(diff_str, diff, vtable, SP_INHOMOGENEOUS, 1);

  // Register the gradient function if we have it.
  if (sp_func_has_deriv(surface1, 1) && sp_func_has_deriv(surface2, 1))
  {
    char diff_grad_str[1024];
    sprintf(diff_grad_str, "Difference gradient"); // FIXME: Yadda
    sp_vtable vtable_g = {.eval = diff_eval_gradient}; // Notice no dtor.
    sp_func_t* diff_grad = sp_func_new(diff_grad_str, diff, vtable_g, SP_INHOMOGENEOUS, 3);
    sp_func_register_deriv(difference, 1, diff_grad);
  }

  return difference;
}

