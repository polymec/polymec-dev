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

#ifndef POLYMEC_POINT_WEIGHT_FUNCTION_H
#define POLYMEC_POINT_WEIGHT_FUNCTION_H

#include "core/polymec.h"
#include "core/point.h"

// A point weight function W(y) is a function that maps the displacement 
// vector y = x - x0 to a scalar value. Usually, it achieves its maximum when 
// x - x0 = 0.
typedef struct point_weight_function_t point_weight_function_t;

// This virtual table must be implemented for any point weight function.
typedef struct
{
  real_t (*value)(void* context, vector_t* y);      // Evaluation of W(y).
  vector_t (*gradient)(void* context, vector_t* y); // Evaluation of grad[W(y)].
  void (*dtor)(void* context);                      // Destructor
} point_weight_function_vtable;

// Creates a point weight function from a name, a context, and a 
// virtual table.
point_weight_function_t* point_weight_function_new(const char* name,
                                                   void* context,
                                                   point_weight_function_vtable vtable);

// Destroys the given point weight function.
void point_weight_function_free(point_weight_function_t* W);

// Returns the name of the point weight function.
const char* point_weight_function_name(point_weight_function_t* W);

// Returns the value of the point weight function for the given displacement 
// vector y = x - x0.
real_t point_weight_function_value(point_weight_function_t* W, vector_t* y);

// Returns the gradient of the point weight function for the given displacement 
// vector y = x - x0.
vector_t point_weight_function_gradient(point_weight_function_t* W, vector_t* y);

#endif
