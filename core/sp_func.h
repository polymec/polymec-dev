#ifndef POLYMEC_SP_FUNC_H
#define POLYMEC_SP_FUNC_H

#include "core/polymec.h"
#include "core/point.h"

#ifdef __cplusplus
extern "C" {
#endif

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
typedef void (*sp_eval_func)(void*, point_t*, double*);

// A destructor for any given context object.
typedef void (*sp_dtor)(void*);

// This virtual table must be implemented by any space-time function.
typedef struct 
{
  sp_eval_func              eval;
  sp_dtor                   dtor;
} sp_vtable;

// Construct a space-time function from the given context, metadata and vtable.
sp_func_t* sp_func_new(const char* name, void* context, sp_vtable vtable,
                       sp_func_homogeneity_t homogeneity,
                       int num_comp);

// Construct a space-time function from a function pointer with the given metadata.
sp_func_t* sp_func_from_func(const char* name, sp_eval_func func, 
                             sp_func_homogeneity_t homogeneity,
                             int num_comp);

// Returns the name of the function.
const char* sp_func_name(sp_func_t* func);

// Returns true if the function is homogeneous in space, false if not.
bool sp_func_is_homogeneous(sp_func_t* func);

// Returns the number of components in the function's range.
int sp_func_num_comp(sp_func_t* func);

// Evaluates the function at the given point, placing the result in result.
void sp_func_eval(sp_func_t* func, point_t* x, double* result);

// Compute the gradient of the function using a centered difference 
// approximation about x0 using the spacing dx.
void sp_func_grad_centered_diff(sp_func_t* func, point_t* x0, vector_t* dx, vector_t* gradient);

// Compute the gradient of the function using Richardson extrapolation 
// applied to a centered difference approximation about x0 using the spacing
// dx.
void sp_func_grad_richardson(sp_func_t* func, point_t* x0, vector_t* dx, vector_t* gradient);

#ifdef __cplusplus
}
#endif

#endif

