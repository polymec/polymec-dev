#ifndef POLYMEC_NONLINEAR_SOLVER_H
#define POLYMEC_NONLINEAR_SOLVER_H

#include "core/polymec.h"
#include "core/adj_graph.h"
#include "core/table.h"

//------------------------------------------------------------------------
//                         nonlinear function
//------------------------------------------------------------------------
// This type represents a nonlinear function F, whose job it is to compute 
// a residual vector R for a given solution vector X. It is a base class that 
// can be subclassed to create actual nonlinear solvers.
typedef struct nonlinear_function_t nonlinear_function_t;

// A function for computing the components of the residual at the given site.
typedef void (*nonlinear_function_eval_residual_func)(void* context, 
                                                      int num_comps, 
                                                      double t, 
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
// residual vector corresponding to those in the solution vector X
// and storing them in R.
void nonlinear_function_eval(nonlinear_function_t* F,
                             double t,
                             double* X,
                             double* R);

//------------------------------------------------------------------------
//                         nonlinear timestepper
//------------------------------------------------------------------------
// This class implements a policy for choosing a timestep given an initial 
// choice, the error history of the nonlinear integration, and any other 
// considerations.
typedef struct nonlinear_timestepper_t nonlinear_timestepper_t;

// A function for computing the size of the timestep, giving a reason for
// the choice. A history of time step sizes (dt), error norms, and iteration 
// counts is stored in the various history arrays, which may be indexed using 
// negative numbers to step backwards in nonlinear iterations. For example, 
// dt_history[0] contains the most recently selected timestep, preceded by 
// dt_history[-1], and so on.
typedef double (*nonlinear_timestepper_compute_dt_func)(void* context, 
                                                        double last_unsuccessful_dt,
                                                        int num_failures,
                                                        double* dt_history,
                                                        double* error_history,
                                                        int* iteration_history,
                                                        int history_length,
                                                        char* explanation);

// A function for determining whether the Jacobian matrix must be recomputed
// during an integration. The histories are available in the same form as 
// described above.
typedef bool (*nonlinear_timestepper_recompute_J_func)(void* context, 
                                                       double* dt_history,
                                                       double* error_history,
                                                       int* iteration_history,
                                                       int history_length);

// A destructor for the timestepper class.
typedef void (*nonlinear_timestepper_dtor)(void* context);

// This virtual table must be implemented by any nonlinear_timestepper_t.
typedef struct 
{
  nonlinear_timestepper_compute_dt_func  compute_dt;
  nonlinear_timestepper_recompute_J_func recompute_J;
  nonlinear_timestepper_dtor             dtor;
} nonlinear_timestepper_vtable;

// Creates a nonlinear timestepper with the given name, context, history 
// length, and virtual table.
nonlinear_timestepper_t* nonlinear_timestepper_new(const char* name, 
                                                   void* context,
                                                   int max_history_length,
                                                   nonlinear_timestepper_vtable vtable);

// Frees a timestepper.
void nonlinear_timestepper_free(nonlinear_timestepper_t* timestepper);

// Returns the (internally-stored) name of the timestepper.
char* nonlinear_timestepper_name(nonlinear_timestepper_t* timestepper);

// Returns the maximum length of the error/iteration history stored in the timestepper.
int nonlinear_timestepper_max_history_length(nonlinear_timestepper_t* timestepper);

// Clears the error/iteration history for the timestepper.
void nonlinear_timestepper_clear_history(nonlinear_timestepper_t* timestepper);

// Returns a new timestep based on the error/iteration history, along with 
// an explanation for the choice.
double nonlinear_timestepper_step_size(nonlinear_timestepper_t* timestepper,
                                       char** explanation);

// Returns true if the Jacobian needs to be recomputed at the present iteration, 
// false if not.
bool nonlinear_timestepper_recompute_jacobian(nonlinear_timestepper_t* timestepper);

// This should be called when a nonlinear iteration converges, with the 
// appropriate time step, error norm, and number of iterations.
void nonlinear_timestepper_converged(nonlinear_timestepper_t* timestepper,
                                     double step_size,
                                     double error_norm,
                                     int num_iterations);

// This should be called when a nonlinear iteration fails to converge with 
// the attempted time step.
void nonlinear_timestepper_failed(nonlinear_timestepper_t* timestepper,
                                  double step_size);

//------------------------------------------------------------------------
//                         nonlinear solver
//------------------------------------------------------------------------
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

// Creates a nonlinear solver that finds the roots of a given nonlinear 
// function at the sites, which depend upon one another according to the 
// given adjacency graph. The solver assumes control over the nonlinear 
// function F and the given timestepper.
nonlinear_solver_t* nonlinear_solver_new(nonlinear_function_t* F,
                                         nonlinear_timestepper_t* timestepper,
                                         adj_graph_t* graph);

// Frees a nonlinear solver.
void nonlinear_solver_free(nonlinear_solver_t* solver);

// Returns the nonlinear function associated with this solver.
nonlinear_function_t* nonlinear_solver_function(nonlinear_solver_t* solver);

// Computes the components of the residual vector R in the nonlinear system 
// at the time t, given the components of the solution vector X.
void nonlinear_solver_compute_residual(nonlinear_solver_t* solver,
                                       double t,
                                       double* X,
                                       double* R);

// Computes the elements of the Jacobian matrix at time t, using the 
// solution vector X and the fractional increment delta. The elements of 
// the Jacobian are stored in the table Jij.
void nonlinear_solver_compute_jacobian(nonlinear_solver_t* solver,
                                       double t,
                                       double* X,
                                       double delta,
                                       double_table_t* Jij);

// Takes a single nonlinear step, returning upon successful convergence, 
// overwriting the values of t and x in place.
void nonlinear_solver_step(nonlinear_solver_t* solver,
                           double* t, double* x);

// Integrates the given nonlinear system from an initial state x1 at time t1 to 
// a final state x2 at t2.
void nonlinear_solver_integrate(nonlinear_solver_t* solver,
                                double t1, double* x1,
                                double t2, double* x2);

// Returns the tolerance for the residual norm that dictates success or 
// failure for a nonlinear iteration. By default, this number is 1e-8.
double nonlinear_solver_tolerance(nonlinear_solver_t* solver);
                                  
// Sets the tolerance for the residual norm that dictates success or 
// failure for a nonlinear iteration.
void nonlinear_solver_set_tolerance(nonlinear_solver_t* solver, double tolerance);

// Returns the absolute maximum number of nonlinear iterations before the 
// solver gives up. By default, this number is 100.
int nonlinear_solver_max_iterations(nonlinear_solver_t* solver);
                                  
// Sets the absolute maximum number of nonlinear iterations before the 
// solver gives up.
void nonlinear_solver_set_max_iterations(nonlinear_solver_t* solver, int max_iters);
                                  
#endif

