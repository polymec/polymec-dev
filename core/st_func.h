// Copyright (c) 2012-2019, Jeffrey N. Johnson
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

/// \addtogroup core core
///@{

/// \class st_func
/// A "space-time" function is an analytic function of space and time.
/// This opaque type encapsulates the notion of such an analytic function
/// and any associated metadata (whether it is homogeneous, constant in
/// time, etc).
/// \refcounted
typedef struct st_func_t st_func_t;

/// \enum st_func_homogeneity_t
/// Enumerated type indicating whether a function is homogeneous in space.
typedef enum
{
  ST_FUNC_HOMOGENEOUS,
  ST_FUNC_HETEROGENEOUS
} st_func_homogeneity_t;

/// \enum st_func_constancy_t
/// Enumerated type indicating whether a function is constant in time.
typedef enum
{
  ST_FUNC_CONSTANT,
  ST_FUNC_NONCONSTANT
} st_func_constancy_t;

/// A function pointer type for evaluating the function at a point.
typedef void (*st_eval_func)(void* context, point_t* x, real_t t, real_t* F);

/// A function pointer type for evaluating the function at n points.
typedef void (*st_eval_n_func)(void* context, point_t* xs, size_t n, real_t t, real_t* Fs);

/// A destructor for any given context object.
typedef void (*st_dtor)(void*);

/// \struct st_func_vtable
/// This virtual table must be implemented by any space-time function.
typedef struct
{
  st_eval_func              eval;
  st_eval_n_func            eval_n;
  st_dtor                   dtor;
} st_func_vtable;

/// Constructs a space-time function from the given context, metadata and vtable.
/// \memberof st_func
st_func_t* st_func_new(const char* name, void* context, st_func_vtable vtable,
                       st_func_homogeneity_t homogeneity,
                       st_func_constancy_t constancy,
                       int num_comp);

/// Constructs a space-time function from a function pointer with the given metadata.
/// The function will be passed NULL as its context.
/// \memberof st_func
st_func_t* st_func_from_func(const char* name, st_eval_func func,
                             st_func_homogeneity_t homogeneity,
                             st_func_constancy_t constancy,
                             int num_comp);

/// Constructs a space_time function from the given spatial function.
/// \memberof st_func
st_func_t* st_func_from_sp_func(sp_func_t* func);

/// Returns the name of the function.
/// \memberof st_func
const char* st_func_name(st_func_t* func);

/// Renames the function. This can be useful for some specialized interfaces.
/// \memberof st_func
void st_func_rename(st_func_t* func, const char* new_name);

/// Returns true if the function is homogeneous in space, false if not.
/// \memberof st_func
bool st_func_is_homogeneous(st_func_t* func);

/// Returns true if the function is constant in time, false if not.
/// \memberof st_func
bool st_func_is_constant(st_func_t* func);

// Returns the number of components in the function's range.
/// \memberof st_func
int st_func_num_comp(st_func_t* func);

/// Returns the context pointer for the given object. Sometimes useful
/// for implementing specialized interfaces.
/// \memberof st_func
void* st_func_context(st_func_t* func);

/// Evaluates the function at the given point, placing the result in result.
/// \memberof st_func
void st_func_eval(st_func_t* func, point_t* x, real_t t, real_t* result);

/// Evaluates the function at n points xs = [x1, x2, ..., xn], placing the
/// result in results = [F1, F2, ..., Fn].
/// \memberof st_func
void st_func_eval_n(st_func_t* func, point_t* xs, size_t n, real_t t, real_t* results);

/// Creates an sp_func from this st_func by "freezing" it at the given time.
/// \memberof st_func
sp_func_t* st_func_freeze(st_func_t* func, real_t t);

/// Constructs a multi-component space-time function from a set of
/// single-valued space-time functions.
/// \memberof st_func
st_func_t* multicomp_st_func_from_funcs(const char* name,
                                        st_func_t** functions,
                                        int num_comp);

/// Constructs a single-component space-time function from one of the
/// components of a multicomponent function.
/// \memberof st_func
st_func_t* st_func_from_component(st_func_t* multicomp_func,
                                  int component);

/// Creates a function that is constant in space and time, with the given
/// components.
/// \memberof st_func
st_func_t* constant_st_func_new(real_t components[], int num_components);

///@}

#endif

