// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_SDT_FUNC_H
#define POLYMEC_SDT_FUNC_H

#include "core/st_func.h"

/// \addtogroup geometry geometry
///@{

/// \class sdt_func
/// This class represents a time-dependent signed distance function, a real-valued function
/// whose zero level set represents a moving surface or shape in space. It is represented by
/// space-time function (st_func) objects.
/// \refcounted
typedef struct sdt_func_t sdt_func_t;

/// \struct sdt_func_vtable
/// This virtual table must be implemented by any time-dependent signed
/// distance function.
typedef struct
{
  real_t (*value)(void* context, point_t* x, real_t t);
  void (*eval_n)(void* context, point_t* xs, size_t n, real_t t, real_t* vals);
  void (*eval_grad)(void* context, point_t* x, real_t t, vector_t* grad);
  void (*eval_n_grad)(void* context, point_t* xs, size_t n, real_t t, vector_t* grads);
  void (*dtor)(void* context);
} sdt_func_vtable;

/// Construct a time-dependent signed distance function with the given name,
/// context, and virtual table.
/// \memberof sdt_func
sdt_func_t* sdt_func_new(const char* name, void* context, sdt_func_vtable vtable);

/// Construct a time-dependent signed distance function with the given name,
/// using the given st_func objects for the distance function and its gradient.
/// \memberof sdt_func
sdt_func_t* sdt_func_from_st_funcs(const char* name, st_func_t* distance, st_func_t* gradient);

/// Returns the name of the signed distance function.
/// \memberof sdt_func
const char* sdt_func_name(sdt_func_t* func);

/// Renames the function. This can be useful for some specialized interfaces.
/// \memberof sdt_func
void sdt_func_rename(sdt_func_t* func, const char* new_name);

/// Returns the context within the sdt_func.
/// \memberof sdt_func
void* sdt_func_context(sdt_func_t* func);

/// Returns the value of the function at the given point x at time t.
/// \memberof sdt_func
real_t sdt_func_value(sdt_func_t* func, point_t* x, real_t t);

/// Evaluates the values of the function at the given n point xs at time t,
/// placing them into vals.
/// \memberof sdt_func
void sdt_func_eval_n(sdt_func_t* func, point_t* xs, size_t n, real_t t, real_t* vals);

/// Computes the gradient of the function at the given point x at time t,
/// placing the result in grad.
/// \memberof sdt_func
void sdt_func_eval_grad(sdt_func_t* func, point_t* x, real_t t, vector_t* grad);

/// Computes the gradients of the function at the given n points xs at time t,
/// placing them in grads.
/// \memberof sdt_func
void sdt_func_eval_n_grad(sdt_func_t* func, point_t* xs, size_t n, real_t t, vector_t* grads);

/// Computes the projection of the point x onto the zero level set of the
/// function at time t, placing the project into proj_x.
/// \memberof sdt_func
void sdt_func_project(sdt_func_t* func, point_t* x, real_t t, point_t* proj_x);

/// Computes the projection of the n points xs onto the zero level set of the
/// function at time t, placing the results into proj_xs.
/// \memberof sdt_func
void sdt_func_project_n(sdt_func_t* func, point_t* xs, size_t n, real_t t, point_t* proj_xs);

///@}


#endif

