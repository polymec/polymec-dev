// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_SDT_FUNC_H
#define POLYMEC_SDT_FUNC_H

#include "core/st_func.h"

// This class represents a time-dependent signed distance function, a real-valued function 
// whose zero level set represents a moving surface or shape in space. It is represented by 
// space-time function (st_func) objects. sdt_func objects are garbage-collected.
typedef struct sdt_func_t sdt_func_t;

// Construct a signed distance function with the given name, using the given 
// st_func objects for the distance function and its gradient.
sdt_func_t* sdt_func_new(const char* name, st_func_t* distance, st_func_t* gradient);

// Returns the name of the signed distance function.
const char* sdt_func_name(sdt_func_t* func);

// Renames the function. This can be useful for some specialized interfaces.
void sdt_func_rename(sdt_func_t* func, const char* new_name);

// Returns the value of the function at the given point x at time t.
real_t sdt_func_eval(sdt_func_t* func, point_t* x, real_t t);

// Computes the gradient of the function at the given point x at time t, 
// placing the result in grad.
void sdt_func_eval_grad(sdt_func_t* func, point_t* x, real_t t, vector_t* grad);

// Computes the projection of the point x onto the zero level set of the 
// function at time t, placing the project into proj_x.
void sdt_func_project(sdt_func_t* func, point_t* x, real_t t, point_t* proj_x);

#endif

