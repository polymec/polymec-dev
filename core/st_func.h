// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_ST_FUNC_H
#define POLYMEC_ST_FUNC_H

#include "core/polymec.h"
#include "core/point.h"
#include "core/sp_func.h"

// A "space-time" function is an analytic function of space and time.
// This opaque type encapsulates the notion of such an analytic function 
// and any associated metadata (whether it is homogeneous, constant in 
// time, etc).
// st_func objects are garbage-collected.
typedef struct st_func_t st_func_t;

// Enumerated type indicating whether a function is homogeneous in space.
typedef enum
{
  ST_FUNC_HOMOGENEOUS,
  ST_FUNC_HETEROGENEOUS
} st_func_homogeneity_t;

// Enumerated type indicating whether a function is constant in time.
typedef enum
{
  ST_FUNC_CONSTANT,
  ST_FUNC_NONCONSTANT
} st_func_constancy_t;

// A function pointer type for evaluating the function at a point.
typedef void (*st_eval_func)(void* context, point_t* x, real_t t, real_t* F);

// A destructor for any given context object.
typedef void (*st_dtor)(void*);

// This virtual table must be implemented by any space-time function.
typedef struct 
{
  st_eval_func              eval;
  st_dtor                   dtor;
} st_func_vtable;

// Constructs a space-time function from the given context, metadata and vtable.
st_func_t* st_func_new(const char* name, void* context, st_func_vtable vtable,
                       st_func_homogeneity_t homogeneity,
                       st_func_constancy_t constancy,
                       int num_comp);

// Constructs a space-time function from a function pointer with the given metadata.
// The function will be passed NULL as its context.
st_func_t* st_func_from_func(const char* name, st_eval_func func, 
                             st_func_homogeneity_t homogeneity,
                             st_func_constancy_t constancy,
                             int num_comp);

// Constructs a space_time function from the given spatial function.
st_func_t* st_func_from_sp_func(sp_func_t* func);

// Returns the name of the function.
const char* st_func_name(st_func_t* func);

// Renames the function. This can be useful for some specialized interfaces.
void st_func_rename(st_func_t* func, const char* new_name);

// Returns true if the function is homogeneous in space, false if not.
bool st_func_is_homogeneous(st_func_t* func);

// Returns true if the function is constant in time, false if not.
bool st_func_is_constant(st_func_t* func);

// Returns the number of components in the function's range.
int st_func_num_comp(st_func_t* func);

// Returns the context pointer for the given object. Sometimes useful 
// for implementing specialized interfaces.
void* st_func_context(st_func_t* func);

// Evaluates the function at the given point, placing the result in result.
void st_func_eval(st_func_t* func, point_t* x, real_t t, real_t* result);

// Registers another function as the nth derivative of this function.
// NOTE: These are vector derivatives, so the first derivative is the 
// NOTE: (3-component) gradient, the second is the (9-component) Hessian,
// NOTE: and so forth.
void st_func_register_deriv(st_func_t* func, int n, st_func_t* nth_deriv);

// Returns true if the nth derivative of this function can be computed (if 
// it was registered), false otherwise.
bool st_func_has_deriv(st_func_t* func, int n);

// Evaluates the nth derivative of this function, placing the result in result.
void st_func_eval_deriv(st_func_t* func, int n, point_t* x, real_t t, real_t* result);

// Returns an st_func representing the nth derivative of this function if 
// available, NULL if not.
st_func_t* st_func_deriv(st_func_t* func, int n);

// Creates an sp_func from this st_func by "freezing" it at the given time.
sp_func_t* st_func_freeze(st_func_t* func, real_t t);

// Constructs a multi-component space-time function from a set of 
// single-valued space-time functions.
st_func_t* multicomp_st_func_from_funcs(const char* name, 
                                        st_func_t** functions,
                                        int num_comp);

// Constructs a single-component space-time function from one of the 
// components of a multicomponent function.
st_func_t* st_func_from_component(st_func_t* multicomp_func,
                                  int component);

// Creates a function that is constant in space and time, with the given 
// components.
st_func_t* constant_st_func_new(real_t components[], int num_components);

#endif

