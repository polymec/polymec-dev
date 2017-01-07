// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

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

// Creates a compactly supported Gaussian weight function with a 
// "shape parameter" epsilon determining the variance of the curve.
point_weight_function_t* gaussian_point_weight_function_new(real_t epsilon);

#endif
