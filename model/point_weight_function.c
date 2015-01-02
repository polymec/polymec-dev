// Copyright (c) 2012-2015, Jeffrey N. Johnson
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

#include "model/point_weight_function.h"

struct point_weight_function_t 
{
  char* name;
  void* context;
  point_weight_function_vtable vtable;
};

point_weight_function_t* point_weight_function_new(const char* name,
                                                   void* context,
                                                   point_weight_function_vtable vtable)
{
  ASSERT(vtable.value != NULL);
  ASSERT(vtable.gradient != NULL);
  point_weight_function_t* W = polymec_malloc(sizeof(point_weight_function_t));
  W->name = string_dup(name);
  W->context = context;
  W->vtable = vtable;
  return W;
}

void point_weight_function_free(point_weight_function_t* W)
{
  if ((W->vtable.dtor != NULL) && (W->context != NULL))
    W->vtable.dtor(W->context);
  string_free(W->name);
  polymec_free(W);
}

const char* point_weight_function_name(point_weight_function_t* W)
{
  return (const char*)W->name;
}

real_t point_weight_function_value(point_weight_function_t* W, vector_t* y)
{
  return W->vtable.value(W->context, y);
}

// Returns the gradient of the point weight function for the given displacement 
// vector y = x - x0.
vector_t point_weight_function_gradient(point_weight_function_t* W, vector_t* y)
{
  return W->vtable.gradient(W->context, y);
}

