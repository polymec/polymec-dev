#ifndef ARBI_st_func_H
#define ARBI_st_func_H

#include "arbi.h"
#include "point.h"
#include "sp_func.h"

#ifdef __cplusplus
extern "C" {
#endif

// A "space-time" function is an analytic function of space and time.
// This opaque type encapsulates the notion of such an analytic function 
// and any associated metadata (whether it is homogeneous, constant in 
// time, etc).
// st_func objects are garbage-collected.
typedef struct st_func_t st_func_t;

// Enumerated type indicating whether a function is homogeneous in space.
typedef enum
{
  ST_HOMOGENEOUS,
  ST_INHOMOGENEOUS
} st_func_homogeneity_t;

// Enumerated type indicating whether a function is constant in time.
typedef enum
{
  ST_CONSTANT,
  ST_NONCONSTANT
} st_func_constancy_t;

// A function pointer type for evaluating the function at a point.
typedef void (*st_eval_func)(void*, point_t*, double, double*);

// A destructor for any given context object.
typedef void (*st_dtor)(void*);

// This virtual table must be implemented by any space-time function.
typedef struct 
{
  st_eval_func              eval;
  st_dtor                   dtor;
} st_vtable;

// Construct a space-time function from the given context, metadata and vtable.
st_func_t* st_func_new(const char* name, void* context, st_vtable vtable,
                       st_func_homogeneity_t homogeneity,
                       st_func_constancy_t constancy,
                       int num_comp);

// Construct a space-time function from a function pointer with the given metadata.
st_func_t* st_func_from_func(const char* name, st_eval_func func, 
                             st_func_homogeneity_t homogeneity,
                             st_func_constancy_t constancy,
                             int num_comp);

// Returns the name of the function.
const char* st_func_name(st_func_t* func);

// Returns true if the function is homogeneous in space, false if not.
bool st_func_is_homogeneous(st_func_t* func);

// Returns true if the function is constant in time, false if not.
bool st_func_is_constant(st_func_t* func);

// Returns the number of components in the function's range.
int st_func_num_comp(st_func_t* func);

// Evaluates the function at the given point, placing the result in result.
void st_func_eval(st_func_t* func, point_t* x, double t, double* result);

// Creates an sp_func from this st_func by "freezing" it at the given time.
sp_func_t* st_func_freeze(st_func_t* func, double t);

#ifdef __cplusplus
}
#endif

#endif

