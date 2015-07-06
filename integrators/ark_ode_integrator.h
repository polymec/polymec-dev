// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_ARK_ODE_INTEGRATOR_H
#define POLYMEC_ARK_ODE_INTEGRATOR_H

#include "integrators/ode_integrator.h"
#include "integrators/newton_pc.h"

// This type of ODE integrator integrates a multi-scale set of ordinary 
// differential equations using semi-implicit Additive Runge Kutte (ARK) 
// methods. The equations integrated take the form
//
// dx/dt = fe(t, x) + fi(t, x)
//
// where fe is a slowly-varying function that can be integrated explicitly, 
// and fi is a quickly-varying function that should be integrated implicitly.

// This type identifies a predictor method to use for implicit integration.
typedef enum
{
  ARK_TRIVIAL_PREDICTOR, // "trivial predictor" -- uses previous time solution
  ARK_MAXORDER_PREDICTOR, // predictor with maximum order polynomial interpolant
  ARK_VARORDER_PREDICTOR, // predictor with variable order polynomial interpolant
  ARK_CUTOFF_PREDICTOR, // predictor with "cutoff" order -- max order in first half of dt,
                        //                                  first order in last half
  ARK_BOOTSTRAP_PREDICTOR // predictor that uses only current step information
} ark_predictor_t;

// This type identifies a Krylov method to use for the Jacobian-Free
// Newton-Krylov version of this integrator.
typedef enum 
{
  JFNK_ARK_GMRES,    // Generalized minimum residual Krylov solver
  JFNK_ARK_FGMRES,   // Flexible GMRES Krylov solver
  JFNK_ARK_BICGSTAB, // Stabilized Biconjugate Gradient Krylov solver
  JFNK_ARK_TFQMR,    // Transpose-Free QMR Krylov solver
  JFNK_ARK_PCG       // Preconditioned Conjugate Gradient Krylov solver
} jfnk_ark_krylov_t;

// This constructs an integrator for the slowly-varying system
// dx/dt = fe(t, x) with the desired order of time accuracy.
// The optional stable_dt_func argument supplies a function that computes the 
// next timestep subject to stability constraints.
ode_integrator_t* explicit_ark_ode_integrator_new(int order, 
                                                  MPI_Comm comm,
                                                  int num_local_values, 
                                                  int num_remote_values, 
                                                  void* context, 
                                                  int (*fe_func)(void* context, real_t t, real_t* x, real_t* fe),
                                                  real_t (*stable_dt_func)(void* context, real_t, real_t* x),
                                                  void (*dtor)(void* context));

// This constructs an integrator for the system with given fe, fi,
// with the desired order of time accuracy, using a fixed-point (functional) 
// iteration method that does not require a Newton-Krylov solver. At least 
// one of fe and fi must be non-NULL. If fi is NULL, the system is assumed to 
// vary "slowly" and will be explicitly integrated; if fe is NULL, the system 
// is assumed to be stiff, and will be implicitly integrated; if both are
// non-NULL, they will be integrated using an adaptive IMEX method.
// max_anderson_accel_dim is the maximum dimension of the underlying Anderson 
// acceleration subspace.
ode_integrator_t* functional_ark_ode_integrator_new(int order, 
                                                    MPI_Comm comm,
                                                    int num_local_values, 
                                                    int num_remote_values, 
                                                    void* context, 
                                                    int (*fe_func)(void* context, real_t t, real_t* x, real_t* fe),
                                                    int (*fi_func)(void* context, real_t t, real_t* x, real_t* fi),
                                                    real_t (*stable_dt_func)(void* context, real_t, real_t* x),
                                                    void (*dtor)(void* context),
                                                    int max_anderson_accel_dim);

// This constructs an integrator for the system with given fe, fi,
// with the desired order of time accuracy, using the Jacobian-Free 
// Newton-Krylov solver of the given type. If the Jacobian-vector 
// product Jy is not given for this integrator, this product will be 
// approximated using difference quotients. At least one of fe and fi must 
// be non-NULL. If fi is NULL, the system is assumed to vary "slowly" and 
// will be explicitly integrated; if fe is NULL, the system is assumed to be 
// stiff, and will be implicitly integrated; if both are non-NULL, they will 
// be integrated using an adaptive IMEX method.
ode_integrator_t* jfnk_ark_ode_integrator_new(int order, 
                                              MPI_Comm comm,
                                              int num_local_values, 
                                              int num_remote_values, 
                                              void* context, 
                                              int (*fe_func)(void* context, real_t t, real_t* x, real_t* fe),
                                              int (*fi_func)(void* context, real_t t, real_t* x, real_t* fi),
                                              real_t (*stable_dt_func)(void* context, real_t, real_t* x),
                                              int (*Jy_func)(void* context, real_t t, real_t* x, real_t* fi, real_t* y, real_t* temp, real_t* Jy),
                                              void (*dtor)(void* context),
                                              newton_pc_t* precond,
                                              jfnk_ark_krylov_t solver_type,
                                              int max_krylov_dim);

// This returns the context pointer passed to the ark_ode_integrator 
// constructor. In general, this will NOT return the same pointer as 
// ode_integrator_context!
void* ark_ode_integrator_context(ode_integrator_t* integrator);

// Sets the (positive) safety factor that is applied to the accuracy-based 
// time step.
void ark_ode_integrator_set_safety_factor(ode_integrator_t* integrator,
                                          real_t factor);

// Sets the method to use for predicting implicit terms in the RK integration.
void ark_ode_integrator_set_predictor(ode_integrator_t* integrator, 
                                      ark_predictor_t predictor);

// Sets the relative and absolute tolerances for integrated quantities.
void ark_ode_integrator_set_tolerances(ode_integrator_t* integrator,
                                       real_t relative_tol, real_t absolute_tol);

// Sets the error weights for evaluating the WRMS norm that is used 
// as a proxy for the quality of the solution. This may be used in lieu of 
// relative and absolute tolerances. Weights are copied into the integrator.
void ark_ode_integrator_set_error_weights(ode_integrator_t* integrator, real_t* weights);

// Sets the error weight function for evaluating the WRMS norm that is used 
// as a proxy for the quality of the solution. This may be used in lieu of 
// relative and absolute tolerances.
void ark_ode_integrator_set_error_weight_function(ode_integrator_t* integrator,
                                                  void (*compute_weights)(void* context, real_t* y, real_t* weights));

// Sets the maximum number of error test failures permitted in attempting 
// a single time step. By default, this value is 7.
void ark_ode_integrator_set_max_err_test_failures(ode_integrator_t* integrator,
                                                  int max_failures);

// Sets the maximum number of nonlinear solver iterations per time step.
// By default, this value is 3.
void ark_ode_integrator_set_max_nonlinear_iterations(ode_integrator_t* integrator,
                                                     int max_iterations);

// Sets the safety factor (coefficient) used in the nonlinear convergence test.
// By default, this value is 0.1.
void ark_ode_integrator_set_nonlinear_convergence_coeff(ode_integrator_t* integrator,
                                                        real_t coefficient);

// Evaluates the explicit part of the system at the given time and with the 
// given solution X, placing the results in fe.
void ark_ode_integrator_eval_fe(ode_integrator_t* integ, real_t t, real_t* X, real_t* fe);

// Evaluates the implicit part of the system at the given time and with the 
// given solution X, placing the results in fi.
void ark_ode_integrator_eval_fi(ode_integrator_t* integ, real_t t, real_t* X, real_t* fi);

// Returns an internal pointer to the preconditioner passed to this 
// integrator during construction time.
newton_pc_t* ark_ode_integrator_preconditioner(ode_integrator_t* integrator);

// Diagnostics for the time integrator.
typedef struct
{
  char* status_message; // borrowed pointer from integrator: do not free.
  long int num_steps, num_explicit_steps, num_accuracy_steps,
           num_step_attempts;
  real_t initial_step_size, last_step_size, next_step_size;
  real_t t;
  long int num_fe_evaluations, num_fi_evaluations;
  long int num_linear_solve_setups;
  long int num_linear_solve_iterations;
  long int num_linear_solve_convergence_failures;
  long int num_error_test_failures;
  long int num_nonlinear_solve_iterations;
  long int num_nonlinear_solve_convergence_failures;
  long int num_preconditioner_evaluations;
  long int num_preconditioner_solves;
} ark_ode_integrator_diagnostics_t;

// Retrieve diagnostics for the time integrator.
void ark_ode_integrator_get_diagnostics(ode_integrator_t* integrator, 
                                        ark_ode_integrator_diagnostics_t* diagnostics);

// Writes time integrator diagnostics to the given file.
void ark_ode_integrator_diagnostics_fprintf(ark_ode_integrator_diagnostics_t* diagnostics, 
                                            FILE* stream);

// This observer type can be used to define objects that respond to actions
// taken by the ark_ode_integrator.
typedef struct ark_ode_observer_t ark_ode_observer_t;

// Creates and returns a newly-allocated observer that observes an 
// am_ode_integrator. Responses are given by the following arguments:
// rhs_computed - This function is called when the right hand side of the ODE
//                system is computed by the integrator.
// Jy_computed - This function is called when the Jacobian-vector product J*y 
//               is computed by the integrator.
// Both of these functions are fed the given context object.
ark_ode_observer_t* ark_ode_observer_new(void* context,
                                         void (*fe_computed)(void* context, real_t t, real_t* x, real_t* rhs),
                                         void (*fi_computed)(void* context, real_t t, real_t* x, real_t* rhs),
                                         void (*Jy_computed)(void* context, real_t t, real_t* x, real_t* rhs, real_t* y, real_t* Jy),
                                         void (*dtor)(void* context));

// Adds the given observer to the given am_ode_integrator. The observer 
// is consumed by the integrator, so no destructor is needed.
void ark_ode_integrator_add_observer(ode_integrator_t* integrator,
                                     ark_ode_observer_t* observer);

#endif

