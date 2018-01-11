// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_SD_FUNC_H
#define POLYMEC_SD_FUNC_H

#include "core/sp_func.h"

// This class represents a signed distance function, a real-valued function whose zero 
// level set represents a surface or shape in space. It is represented by spatial function 
// (sp_func) objects. sd_func objects are garbage-collected.
typedef struct sd_func_t sd_func_t;

// This virtual table must be implemented by any signed distance function.
typedef struct 
{
  real_t (*value)(void* context, point_t* x);
  void (*eval_n)(void* context, point_t* xs, size_t n, real_t* vals);
  void (*eval_grad)(void* context, point_t* x, vector_t* grad);
  void (*eval_n_grad)(void* context, point_t* xs, size_t n, vector_t* grads);
  void (*dtor)(void* context);
} sd_func_vtable;

// Construct a signed distance function with the given name, context, and 
// virtual table.
sd_func_t* sd_func_new(const char* name, void* context, sd_func_vtable vtable);

// Construct a signed distance function with the given name, using the given 
// sp_func objects for the distance function and its gradient.
sd_func_t* sd_func_from_sp_funcs(const char* name, sp_func_t* distance, sp_func_t* gradient);

// Returns the name of the signed distance function.
const char* sd_func_name(sd_func_t* func);

// Renames the function. This can be useful for some specialized interfaces.
void sd_func_rename(sd_func_t* func, const char* new_name);

// Returns the context within the sd_func.
void* sd_func_context(sd_func_t* func);

// Returns the value of the function at the given point x.
real_t sd_func_value(sd_func_t* func, point_t* x);

// Evaluates the values of the function at the given n point xs, placing 
// them into vals.
void sd_func_eval_n(sd_func_t* func, point_t* xs, size_t n, real_t* vals);

// Computes the gradient of the function at the given point x, placing the 
// result in grad.
void sd_func_eval_grad(sd_func_t* func, point_t* x, vector_t* grad);

// Computes the gradients of the function at the given n points xs, placing 
// them in grads.
void sd_func_eval_n_grad(sd_func_t* func, point_t* xs, size_t n, vector_t* grads);

// Computes the projection of the point x onto the zero level set of the 
// function, placing the result into proj_x.
void sd_func_project(sd_func_t* func, point_t* x, point_t* proj_x);

// Computes the projection of the n points xs onto the zero level set of the 
// function, placing the results into proj_xs.
void sd_func_project_n(sd_func_t* func, point_t* xs, size_t n, point_t* proj_xs);

#endif

