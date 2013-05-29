#ifndef POLYMEC_NONLINEAR_SOLVER_H
#define POLYMEC_NONLINEAR_SOLVER_H

#include "core/polymec.h"
#include "core/table.h"

// This class solves a nonlinear system F(X) = R, where F is a nonlinear 
// function mapping a solution vector X to a residual vector R. The system 
// itself has a specified dimension that is the size of the vectors X and R.
// In the context of partial differential equations, it is sometimes 
// advantageous to specify a number of "components" to a solution vector X 
// that describes the number of unknowns corresponding to a single "site" at
// which the nonlinear function may be evaluated. In this parlance, the 
// dimension is the product of the number of components and the number of 
// sites.
typedef struct nonlinear_solver_t nonlinear_solver_t;

// This type represents a nonlinear function F, whose job it is to compute 
// a residual vector R for a given solution vector X. It is a base class that 
// can be subclassed to create actual nonlinear solvers.
typedef struct nonlinear_function_t nonlinear_function_t;

// A function for computing the components of the residual at the given site.
typedef void (*nonlinear_function_eval_residual_func)(void* context, 
                                                      int num_comps, 
                                                      int site, 
                                                      double* x, 
                                                      double* R);

// A destructor for nonlinear_functions.
typedef void (*nonlinear_function_dtor)(void* context);

// This virtual table must be implemented by any nonlinear_function_t.
typedef struct 
{
  nonlinear_function_eval_residual_func   eval_residual;
  nonlinear_function_dtor                 dtor;
} nonlinear_function_vtable;

// Creates a nonlinear function with the given name, context, 
// number of solution components, and virtual table.
nonlinear_function_t* nonlinear_function_new(const char* name, 
                                             void* context,
                                             int num_comps,
                                             nonlinear_function_vtable vtable);

// Frees a nonlinear function.
void nonlinear_function_free(nonlinear_function_t* F);

// Returns the name of the nonlinear function (internally stored).
char* nonlinear_function_name(nonlinear_function_t* F);

// Returns the context object for this nonlinear function.
void* nonlinear_function_context(nonlinear_function_t* F);

// Returns the number of components per site in this nonlinear function.
int nonlinear_function_num_comps(nonlinear_function_t* F);

// Evaluates the nonlinear function, computing the components of the 
// residual vector corresponding to those in the solution vector X for a 
// given site and storing them in R.
void nonlinear_function_eval(nonlinear_function_t* F,
                             int site,
                             double* X,
                             double* R);

// Creates a nonlinear solver that finds the roots of a given nonlinear 
// function at a given number of sites. The solver assumes control over the 
// nonlinear function F.
nonlinear_solver_t* nonlinear_solver_new(nonlinear_function_t* F,
                                         int num_sites);

// Frees a nonlinear solver.
void nonlinear_solver_free(nonlinear_solver_t* solver);

// Returns the nonlinear function associated with this solver.
nonlinear_function_t* nonlinear_solver_function(nonlinear_solver_t* solver);

// Computes the components of the residual vector R at the given site in 
// the nonlinear solver, given the components of the solution vector X.
void nonlinear_solver_compute_residual(nonlinear_solver_t* solver,
                                       int site,
                                       double* X,
                                       double* R);

// Computes the elements of the Jacobian matrix at the site i, using the 
// solution vector X and the fractional increment delta. The elements of 
// the Jacobian are stored in the table Jij.
void nonlinear_solver_compute_jacobian(nonlinear_solver_t* solver,
                                       int site,
                                       double* X,
                                       double delta,
                                       double_table_t* Jij);

#endif

