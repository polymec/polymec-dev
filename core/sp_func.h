// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_SP_FUNC_H
#define POLYMEC_SP_FUNC_H

#include "core/polymec.h"
#include "core/point.h"

/// \addtogroup core core
///@{

/// \class sp_func
/// A "space" function (or "spatial" function) is an analytic function of
/// space only. This opaque type encapsulates the notion of such an analytic
/// function and any associated metadata (whether it is homogeneous or not, etc).
/// \refcounted
typedef struct sp_func_t sp_func_t;

/// \enum sp_func_homogeneity_t
/// Enumerated type indicating whether a function is homogeneous.
typedef enum
{
  SP_FUNC_HOMOGENEOUS,
  SP_FUNC_HETEROGENEOUS
} sp_func_homogeneity_t;

/// A function pointer type for evaluating the function at a point.
typedef void (*sp_eval_func)(void* context, point_t* x, real_t* F);

/// A function pointer type for evaluating the function at n points.
typedef void (*sp_eval_n_func)(void* context, point_t* xs, size_t n, real_t* Fs);

/// A destructor for any given context object.
typedef void (*sp_dtor)(void*);

/// \struct sp_func_vtable
/// This virtual table must be implemented by any space-time function.
typedef struct
{
  sp_eval_func              eval;
  sp_eval_n_func            eval_n;
  sp_dtor                   dtor;
} sp_func_vtable;

/// Constructs a space-time function from the given context, metadata and vtable.
/// \memberof sp_func
sp_func_t* sp_func_new(const char* name, void* context, sp_func_vtable vtable,
                       sp_func_homogeneity_t homogeneity,
                       int num_comp);

/// Constructs a space-time function from a function pointer with the given metadata.
/// The function will be passed NULL as its context.
/// \memberof sp_func
sp_func_t* sp_func_from_func(const char* name, sp_eval_func func,
                             sp_func_homogeneity_t homogeneity,
                             int num_comp);

/// Returns the name of the function.
/// \memberof sp_func
const char* sp_func_name(sp_func_t* func);

/// Renames the function. This can be useful for some specialized interfaces.
/// \memberof sp_func
void sp_func_rename(sp_func_t* func, const char* new_name);

/// Returns true if the function is homogeneous in space, false if not.
/// \memberof sp_func
bool sp_func_is_homogeneous(sp_func_t* func);

/// Returns the number of components in the function's range.
/// \memberof sp_func
int sp_func_num_comp(sp_func_t* func);

/// Returns the context pointer for the given object. Sometimes useful
/// for implementing specialized interfaces.
/// \memberof sp_func
void* sp_func_context(sp_func_t* func);

/// Evaluates the function at the given point, placing the result in result.
/// \memberof sp_func
void sp_func_eval(sp_func_t* func, point_t* x, real_t* result);

/// Evaluates the function at the given n points xs = {x1, x2, ..., xn},
/// placing the result in result = {F1, F2, ..., Fn}.
/// \memberof sp_func
void sp_func_eval_n(sp_func_t* func, point_t* xs, size_t n, real_t* result);

/// Creates a function that is constant everywhere in space, with the given
/// components.
/// \memberof sp_func
sp_func_t* constant_sp_func_new(real_t components[], int num_components);

///@}

#endif

