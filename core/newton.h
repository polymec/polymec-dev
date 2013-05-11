#ifndef POLYMEC_NEWTON_H
#define POLYMEC_NEWTON_H

#include "core/polymec.h"

// This defines a function that computes the value and derivative of a 
// (real-valued) function of a single variable given a context.
typedef void (*nonlinear_function_t)(void*, double, double*, double*);

// Solve the nonlinear equation F(x) = 0 using Newton-Raphson iteration on 
// the interval [min, max] using the initial estimate x, with the specified 
// maximum number of iterations and desired tolerance.
// Returns true if the solution was achieved, false otherwise.
bool newton_solve(nonlinear_function_t F, void* context, double* x, double min, double max, double tolerance, int max_iters); 

// Solve the nonlinear equation F(x) = 0 using Brent's method on 
// the interval [x1, x2] with the specified maximum number of iterations 
// and desired tolerance. Returns the solution x. NOTE that the nonlinear
// function F need not have a derivative for Brent's method.
double brent_solve(nonlinear_function_t F, void* context, double x1, double x2, double tolerance, int max_iters, double* error); 

// This defines a vector-valued function that computes its value 
// given a context. 
typedef void (*nonlinear_vector_function_t)(void*, double*, double*);

// This function prototype describes functions that compute Jacobians.
typedef void (*jacobian_function_t)(void*, nonlinear_vector_function_t, int, double*, double*);

// This struct represents a nonlinear system of equations F(x) = 0.
typedef struct
{
  int dim;                                // Dimension of the system.
  nonlinear_vector_function_t compute_F;  // Vector-valued function F.
  jacobian_function_t compute_J;          // Jacobian calculation (optional).
  void* context;                          // Context pointer.
} nonlinear_system_t;

// Solve the nonlinear system F(x) = 0 using a globally-convergent version 
// of Newton's method with the initial estimate vector x, with the specified 
// maximum number of iterations and desired tolerance. Returns true if the 
// solution converged, false otherwise. If check_solution is set to true, it is possible that spurious convergence occurred.
bool newton_solve_system(nonlinear_system_t* system, double* x, double tolerance, int max_iters, int* num_iters);

// Solve the nonlinear system F(x) = 0 using Broyden's method with the 
// initial estimate vector x, with the specified maximum number of iterations 
// and desired tolerance. Returns true if the solution converged, false otherwise. 
// If check_solution is set to true, it is possible that spurious convergence occurred.
bool broyden_solve_system(nonlinear_system_t* system, double* x, double tolerance, int max_iters, int* num_iters);

#endif

