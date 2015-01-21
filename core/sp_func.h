// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_SP_FUNC_H
#define POLYMEC_SP_FUNC_H

#include "core/polymec.h"
#include "core/point.h"

// A "space" function (or "spatial" function) is an analytic function of 
// space only. This opaque type encapsulates the notion of such an analytic 
// function and any associated metadata (whether it is homogeneous or not, etc).
// sp_func objects are garbage-collected.
typedef struct sp_func_t sp_func_t;

// Enumerated type indicating whether a function is homogeneous.
typedef enum
{
  SP_HOMOGENEOUS,
  SP_INHOMOGENEOUS
} sp_func_homogeneity_t;

// A function pointer type for evaluating the function at a point.
typedef void (*sp_eval_func)(void* context, point_t* x, real_t* F);

// A function pointer type for evaluating the nth derivative of the 
// function at a point.
typedef void (*sp_eval_deriv_func)(void*, int, point_t*, real_t*);

// A function pointer type for returning whether a function has a 
// derivative. This must be supplied if eval_deriv is given.
typedef bool (*sp_has_deriv_func)(void*, int);

// A destructor for any given context object.
typedef void (*sp_dtor)(void*);

// This virtual table must be implemented by any space-time function.
typedef struct 
{
  sp_eval_func              eval;
  sp_eval_deriv_func        eval_deriv; // Optional
  sp_has_deriv_func         has_deriv; // Optional, but must be given if eval_deriv is given
  sp_dtor                   dtor;
} sp_vtable;

// Construct a space-time function from the given context, metadata and vtable.
sp_func_t* sp_func_new(const char* name, void* context, sp_vtable vtable,
                       sp_func_homogeneity_t homogeneity,
                       int num_comp);

// Construct a space-time function from a function pointer with the given metadata.
// The function will be passed NULL as its context.
sp_func_t* sp_func_from_func(const char* name, sp_eval_func func, 
                             sp_func_homogeneity_t homogeneity,
                             int num_comp);

// Returns the name of the function.
const char* sp_func_name(sp_func_t* func);

// Renames the function. This can be useful for some specialized interfaces.
void sp_func_rename(sp_func_t* func, const char* new_name);

// Returns true if the function is homogeneous in space, false if not.
bool sp_func_is_homogeneous(sp_func_t* func);

// Returns the number of components in the function's range.
int sp_func_num_comp(sp_func_t* func);

// Returns the context pointer for the given object. Sometimes useful 
// for implementing specialized interfaces.
void* sp_func_context(sp_func_t* func);

// Evaluates the function at the given point, placing the result in result.
void sp_func_eval(sp_func_t* func, point_t* x, real_t* result);

// Registers another function as the nth derivative of this function.
// NOTE: These are vector derivatives, so the first derivative is the 
// NOTE: (3-component) gradient, the second is the (9-component) Hessian,
// NOTE: and so forth.
void sp_func_register_deriv(sp_func_t* func, int n, sp_func_t* nth_deriv);

// Returns true if the nth derivative of this function can be computed (if 
// it was registered), false otherwise.
bool sp_func_has_deriv(sp_func_t* func, int n);

// Evaluates the nth derivative of this function, placing the result in result.
void sp_func_eval_deriv(sp_func_t* func, int n, point_t* x, real_t* result);

// Creates a function that is constant everywhere in space, with the given 
// components.
sp_func_t* constant_sp_func_new(real_t components[], int num_components);

#endif

