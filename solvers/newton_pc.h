// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_NEWTON_PC_H
#define POLYMEC_NEWTON_PC_H

#include "core/polymec.h"

/// \addtogroup solvers solvers
///@{

/// \class newton_pc
/// This type represents a preconditioner for Newton methods on
/// large distributed systems, used for:
/// 1. Finding the solution x of a nonlinear system F(t, x) = 0.
/// 2. Integrating a nonlinear system of ODEs dx/dt = F(t, x).
/// 3. Solving a system of nonlinear differential algebraic equations (DAE):
///     F(t, x, xdot) = 0.
typedef struct newton_pc_t newton_pc_t;

/// \struct newton_pc_vtable
/// Every Newton preconditioner has to implement the functions in the 
/// following virtual table.
typedef struct
{
  /// Method to reset the preconditioner to its original state at time t.
  void (*reset)(void* context, real_t t);

  /// Method for computing the preconditioner operator 
  /// P = alpha * I + beta * dF/dx + gamma * dF/d(xdot).
  /// where F(t, x, xdot) is a nonlinear function whose solution x is sought.
  void (*compute_p)(void* context, 
                    real_t alpha, real_t beta, real_t gamma, 
                    real_t t, real_t* x, real_t* xdot);

  /// Method to solve the preconditioner equation P * z = r for z, given the 
  /// state (t, x, xdot) and a residual L2 norm tolerance. The error L2 norm 
  /// should be stored in *error_L2_norm.
  bool (*solve)(void* context, 
                real_t t, real_t* x, real_t* xdot, real_t tolerance,
                real_t* r, real_t* z, real_t* error_L2_norm);

  /// Destructor.
  void (*dtor)(void* context);
} newton_pc_vtable;

/// \enum newton_pc_side_t
/// Preconditioner sided-ness.
typedef enum
{
  NEWTON_PC_LEFT,
  NEWTON_PC_RIGHT,
  NEWTON_PC_BOTH
} newton_pc_side_t;

/// Constructs a preconditioner that represents the linear combination
/// alpha * I + beta * dF/dx + gamma * dF/dx, where F(t, x[, xdot]) = 0 is a 
/// function representing a nonlinear ODE [or differential/algebraic) system, 
/// I is the identity matrix, and dF/dx [and dF/dxdot] are the derivatives of 
/// F with respect to x [and xdot]. The last argument, side, indicates whether 
/// the preconditioner should be applied to the left, right, or both sides 
/// of the equation.
/// \memberof newton_pc
newton_pc_t* newton_pc_new(const char* name,
                           void* context,
                           newton_pc_vtable vtable,
                           newton_pc_side_t side);

/// Frees the preconditioner.
/// \memberof newton_pc
void newton_pc_free(newton_pc_t* precond);

/// Returns an internal string containing the name of the preconditioner.
/// \memberof newton_pc
char* newton_pc_name(newton_pc_t* precond);

/// Returns an internal pointer to the context originally given to this 
/// Newton preconditioner. Use at your own risk.
/// \memberof newton_pc
void* newton_pc_context(newton_pc_t* precond);

/// Returns the "sidedness" of the Newton preconditioner.
/// \memberof newton_pc
newton_pc_side_t newton_pc_side(newton_pc_t* precond);

/// Resets the preconditioner to its original state at time t.
/// \memberof newton_pc
void newton_pc_reset(newton_pc_t* precond, real_t t);

/// Sets the tolerance to which the L2 norm of r - P*z will be compared in 
/// order for newton_pc_solve to succeed.
/// \memberof newton_pc
void newton_pc_set_tolerance(newton_pc_t* precond, real_t tolerance);

/// Performs any setup required to computes the preconditioner matrix at the 
/// point (t, x, xdot) in solution space, computing 
/// alpha * I + beta * dF/dx + gamma * dF/d(xdot) 
/// for the given alpha, beta, and gamma. Note that xdot should be NULL if F
/// is only a function of x.
/// \memberof newton_pc
void newton_pc_setup(newton_pc_t* precond, 
                     real_t alpha, real_t beta, real_t gamma,
                     real_t t, real_t* x, real_t* xdot);

/// Solves the preconditioner system P*z = r, placing the solution in z.
/// The time, solution, and (possibly) its time derivative are given.
/// This can only be called after newton_pc_setup has been called.
/// Returns true if the solve succeeded to within the preconditioner's 
/// given tolerance, false otherwise.
/// \memberof newton_pc
bool newton_pc_solve(newton_pc_t* precond, 
                     real_t t, real_t* x, real_t* xdot,
                     real_t* r, real_t* z);

// The following "fix/unfix" methods are for fixing and unfixing the alpha, 
// beta, and gamma coefficients, and should really only be used in tandem with 
// the newton_solver, and not Newton-Krylov ODE/DAE solvers.

/// Fixes the alpha, beta, gamma coefficients to fixed values so that newton_pc_setup always 
/// uses these fixed values instead of those requested.
/// \memberof newton_pc
void newton_pc_fix_coefficients(newton_pc_t* precond, real_t alpha0, real_t beta0, real_t gamma0);

/// Unfixes the alpha, beta, gamma coefficients if they were previously fixed. 
/// Has no effect if these coefficients weren't previously fixed.
/// \memberof newton_pc
void newton_pc_unfix_coefficients(newton_pc_t* precond);

/// Returns true if the alpha, beta, gamma coefficients have been fixed, false 
/// if not.
/// \memberof newton_pc
bool newton_pc_coefficients_fixed(newton_pc_t* precond);

///@}

#endif

