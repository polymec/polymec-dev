#ifndef POLYMEC_NEWTON_H
#define POLYMEC_NEWTON_H

#include "core/polymec.h"

#ifdef __cplusplus
extern "C" {
#endif

// This defines a function that computes the value and derivative of a 
// (real-valued) function of a single variable given a context.
typedef void (*nonlinear_function_t)(void*, double, double*, double*);

// This defines a vector-valued function that computes its value and 
// Jacobian given a context. 
typedef void (*nonlinear_vector_function_t)(void*, double*, double*, double*);

// Solve the nonlinear equation F(x) = 0 using Newton-Raphson iteration on 
// the interval [min, max] using the initial estimate x, with the specified 
// maximum number of iterations and desired tolerance.
// Returns true if the solution was achieved, false otherwise.
bool newton_solve(nonlinear_function_t F, void* context, double* x, double min, double max, double tolerance, int max_iters); 

// Solve the nonlinear equation F(x) = 0 using Picard iteration 
// using the initial estimate x, with the specified maximum number of 
// iterations and desired tolerance.
// Returns true if the solution was achieved, false otherwise.
bool picard_solve(nonlinear_function_t F, void* context, double* x, double tolerance, int max_iters); 

// Solve the nonlinear system F(x) = 0 using Picard iteration using the 
// initial estimate vector x, with the specified maximum number of iterations 
// and desired tolerance.
// Returns true if the solution was achieved, false otherwise.
bool picard_solve_system(int dim, nonlinear_vector_function_t F, void* context, double* x, double tolerance, int max_iters); 

// Solve the nonlinear equation F(x) = 0 using a hybrid Newton-Raphson/Picard 
// iteration on the interval [min, max] using the initial estimate x, with 
// the desired tolerance. The Newton iteration is tried first, and performs 
// newton_iters iterations before switching to a Picard iteration, which 
// performs picard_iters iterations. This cycle is performed num_cycles times.
// Returns true if the solution was achieved, false otherwise.
bool newton_picard_solve(nonlinear_function_t F, void* context, double* x, double min, double max, double tolerance, int newton_iters, int picard_iters, int num_cycles); 

// Solve the nonlinear equation F(x) = 0 using Brent's method on 
// the interval [x1, x2] with the specified maximum number of iterations 
// and desired tolerance. Returns the solution x. NOTE that the nonlinear
// function F need not have a derivative for Brent's method.
double brent_solve(nonlinear_function_t F, void* context, double x1, double x2, double tolerance, int max_iters, double* error); 

#ifdef __cplusplus
}
#endif

#endif

