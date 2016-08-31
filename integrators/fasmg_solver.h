// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_FASMG_SOLVER_H
#define POLYMEC_FASMG_SOLVER_H

#include "core/polymec.h"

// The fasmg_solver class provides an interface for solving nonlinear systems 
// of equations using the Full Approximation Storage (FAS) multigrid method.
// A FAS multigrid method consists of the following algorithms:
//
// 1. A nonlinear operator A(X) (fasmg_operator), defined on a given 
//    discretization and operating on a solution vector X on that discretization.
//    This operator also implements a nonlinear relaxation scheme that is used
//    to smooth the error in A(X) = B.
// 2. A coarsening scheme (fasmg_coarsener) that creates a coarser version of 
//    a given discretization.
// 3. A restriction scheme (fasmg_restrictor) that projects a vector from a 
//    fine discretization to a coarser one.
// 4. A prolongation scheme (fasmg_prolongator) that interpolates a vector 
//    from a coarser discretization to a finer one.
// 5. A cycle scheme (fasmg_cycle) that implements and performs a multigrid 
//    cycle (such as a V cycle, a W cycle, or a full multigrid (FMG) cycle).
//
// The fasmg_solver itself manages these algorithms, using them to solve 
// A(X) = B for a given X on a given discretization. Note that the discretization 
// here is considered a black box, and is referred to only as a void*. The 
// various algorithms can use and manipulate discretizations of this type, but 
// the multigrid solver itself is neutral to the representation of the 
// discretization.
typedef struct fasmg_solver_t fasmg_solver_t;

// This class represents a discretization to be used with the fasmg_solver.
// Usually it is an opaque container for an application-specific data structure.
typedef struct fasmg_grid_t fasmg_grid_t;

// This class represents a nonlinear operator A(X).
typedef struct fasmg_operator_t fasmg_operator_t;

// This class represents a scheme for creating coarsened discretizations
// from an existing one.
typedef struct fasmg_coarsener_t fasmg_coarsener_t;

// This class represents a scheme for restricting a solution vector on a fine 
// discretization to one on a coarser discretization.
typedef struct fasmg_restrictor_t fasmg_restrictor_t;

// This class represents a scheme for prolongating a solution vector on a coarse 
// discretization to one on a finer discretization.
typedef struct fasmg_prolongator_t fasmg_prolongator_t;

// This class represents a multigrid cycle.
typedef struct fasmg_cycle_t fasmg_cycle_t;

// Creates a fasmg_solver that uses the given algorithms to solve A(X) = B.
// The solver assumes responsibility for managing all of the algorithms 
// and freeing their resources.
fasmg_solver_t* fasmg_solver_new(fasmg_operator_t* A,
                                 fasmg_coarsener_t* coarsener,
                                 fasmg_restrictor_t* restrictor,
                                 fasmg_prolongator_t* prolongator,
                                 fasmg_cycle_t* cycle);

// Frees the given fasmg_solver.
void fasmg_solver_free(fasmg_solver_t* solver);

// Sets the maximum number of multigrid cycles that can be used to solve
// a system.
void fasmg_solver_set_max_cycles(fasmg_solver_t* solver,
                                 int max_cycles);

// Sets the maximum residual norm allowed for a solution to be considered 
// convergent.
void fasmg_solver_set_max_residual_norm(fasmg_solver_t* solver,
                                        real_t max_residual_norm);

// Returns a newly-constructed representation of the discretization that can be 
// used with this FASMG solver, created from the given underlying representation,
// with the given number of degrees of freedom, and using the given destructor 
// to free resources for grids created by the coarsener. The given pointer is 
// not managed by the solver, however.
fasmg_grid_t* fasmg_solver_grid(fasmg_solver_t* solver,
                                void* discretization,
                                size_t num_dof,
                                void (*discretization_dtor)(void*));

// Solves the system A(X) = B on the given discretization for the given 
// vectors B and X. Returns true if the solution converged, false if not. The 
// residual L2 norm and the number of cycles are stored in their respective 
// variables.
bool fasmg_solver_solve(fasmg_solver_t* solver,
                        fasmg_grid_t* grid,
                        real_t* B,
                        real_t* X,
                        real_t* residual_norm,
                        int* num_cycles);

// Performs a single multigrid cycle for A(X) = B on the given discretization.
void fasmg_solver_cycle(fasmg_solver_t* solver,
                        fasmg_grid_t* grid,
                        real_t* B,
                        real_t* X);

//------------------------------------------------------------------------
//                            fasmg_grid
//------------------------------------------------------------------------

// fasmg_grid objects are created by the multigrid solver as needed, so there 
// are no explicit grid constructors.

// Destroys this fasmg_grid.
void fasmg_grid_free(fasmg_grid_t* grid);

// Returns the adjacent coarser grid in the multigrid hierarchy, or NULL if 
// there is no such grid.
fasmg_grid_t* fasmg_grid_coarser(fasmg_grid_t* grid);

// Returns the adjacent finer grid in the multigrid hierarchy, or NULL if 
// there is no such grid.
fasmg_grid_t* fasmg_grid_finer(fasmg_grid_t* grid);

//------------------------------------------------------------------------
//                            fasmg_operator
//------------------------------------------------------------------------

// This virtual table determines the implementation of the fasmg_operator.
typedef struct
{
  // Applies the operator to a vector X on a grid, storing A(X) in AX.
  void (*apply)(void* context, void* grid, real_t* X, real_t* AX);
  // Updates X by performing a relaxation step on A(X) = B.
  void (*relax)(void* context, void* grid, real_t* B, real_t* X);
  // Destructor.
  void (*dtor)(void* context);
} fasmg_operator_vtable;

// Creates a fasmg_operator that uses the given context and vtable to 
// implement its behavior. 
fasmg_operator_t* fasmg_operator_new(const char* name, 
                                     void* context,
                                     fasmg_operator_vtable vtable);

// Frees the given fasmg_operator.
void fasmg_operator_free(fasmg_operator_t* A);

// Returns the name of the operator (internally stored).
char* fasmg_operator_name(fasmg_operator_t* A);

// Returns the context object for this integrator.
void* fasmg_operator_context(fasmg_operator_t* A);

// Applies the operator A to the vector X, computing A(X) and storing it in AX.
void fasmg_operator_apply(fasmg_operator_t* A,
                          fasmg_grid_t* grid,
                          real_t* X, 
                          real_t* AX);

// Performs a relaxation step, updating X to eliminate high-frequency components 
// of the error in A(X) = B.
void fasmg_operator_relax(fasmg_operator_t* A,
                          fasmg_grid_t* grid,
                          real_t* B,
                          real_t* X);

// Computes the residual R = B - A(X) for this operator.
void fasmg_operator_compute_residual(fasmg_operator_t* A,
                                     fasmg_grid_t* grid,
                                     real_t* B,
                                     real_t* X,
                                     real_t* R);

//------------------------------------------------------------------------
//                            fasmg_coarsener
//------------------------------------------------------------------------

// This virtual table determines the implementation of the fasmg_coarsener.
typedef struct
{
  // This function returns true if the given grid can be coarsened further, 
  // false if not.
  bool (*can_coarsen)(void* context, void* grid);
  // This function returns a new, coarser grid, given a fine grid. It also 
  // stores the coarse number of degrees of freedom in coarse_dof
  void* (*coarser_grid)(void* context, void* grid, size_t* coarse_dof);
  // Destructor.
  void (*dtor)(void* context);
} fasmg_coarsener_vtable;

// Creates a fasmg_coarsener that uses the given context and vtable to 
// implement its behavior. 
fasmg_coarsener_t* fasmg_coarsener_new(const char* name, 
                                       void* context,
                                       fasmg_coarsener_vtable vtable);

// Frees the given fasmg_coarsener.
void fasmg_coarsener_free(fasmg_coarsener_t* coarsener);

// Returns the name of the coarsener (internally stored).
char* fasmg_coarsener_name(fasmg_coarsener_t* coarsener);

// Returns the context object for this integrator.
void* fasmg_coarsener_context(fasmg_coarsener_t* coarsener);

// Returns true if the coarsener can coarsen the given grid, false if not.
bool fasmg_coarsener_can_coarsen(fasmg_coarsener_t* coarsener,
                                 fasmg_grid_t* grid);

// Creates a series of increasingly coarser grids associated with the given
// one, terminating when fasmg_coarsener_can_coarsen(grid) returns false.
// The coarser grids are accessible via fasmg_grid_coarser(grid).
void fasmg_coarsener_coarsen(fasmg_coarsener_t* coarsener,
                             fasmg_grid_t* grid);

//------------------------------------------------------------------------
//                            fasmg_restrictor
//------------------------------------------------------------------------

// This virtual table determines the implementation of the fasmg_restrictor.
typedef struct
{
  // Projects the given vector on a fine grid to a new vector on a coarse 
  // grid.
  void (*project)(void* context, void* fine_grid, real_t* fine_X, void* coarse_grid, real_t* coarse_X);
  // Destructor.
  void (*dtor)(void* context);
} fasmg_restrictor_vtable;

// Creates a fasmg_restrictor that uses the given context and vtable to 
// implement its behavior. 
fasmg_restrictor_t* fasmg_restrictor_new(const char* name, 
                                         void* context,
                                         fasmg_restrictor_vtable vtable);

// Frees the given fasmg_restrictor.
void fasmg_restrictor_free(fasmg_restrictor_t* restrictor);

// Returns the name of the restrictor (internally stored).
char* fasmg_restrictor_name(fasmg_restrictor_t* restrictor);

// Returns the context object for this integrator.
void* fasmg_restrictor_context(fasmg_restrictor_t* restrictor);

// Projects the vector fine_X, defined on fine_grid, to the 
// vector coarse_X, defined on the coarser grid corresponding to 
// fine_grid.
void fasmg_restrictor_project(fasmg_restrictor_t* restrictor,
                              fasmg_grid_t* fine_grid,
                              real_t* fine_X,
                              real_t* coarse_X);

//------------------------------------------------------------------------
//                            fasmg_prolongator
//------------------------------------------------------------------------

// This virtual table determines the implementation of the fasmg_prolongator.
typedef struct
{
  // Interpolates the given vector on a fine grid to a new vector on a coarse 
  // grid.
  void (*interpolate)(void* context, void* coarse_grid, real_t* coarse_X, void* fine_grid, real_t* fine_X);
  // Destructor.
  void (*dtor)(void* context);
} fasmg_prolongator_vtable;

// Creates a fasmg_prolongator that uses the given context and vtable to 
// implement its behavior. 
fasmg_prolongator_t* fasmg_prolongator_new(const char* name, 
                                           void* context,
                                           fasmg_prolongator_vtable vtable);

// Frees the given fasmg_prolongator.
void fasmg_prolongator_free(fasmg_prolongator_t* prolongator);

// Returns the name of the prolongator (internally stored).
char* fasmg_prolongator_name(fasmg_prolongator_t* prolongator);

// Returns the context object for this integrator.
void* fasmg_prolongator_context(fasmg_prolongator_t* prolongator);

// Interpolates the vector coarse_X, defined on coarse_grid, to the 
// vector fine_X, defined on the finer grid corresponding to coarse_grid.
void fasmg_prolongator_interpolate(fasmg_prolongator_t* prolongator,
                                   fasmg_grid_t* coarse_grid,
                                   real_t* coarse_X,
                                   real_t* fine_X);

//------------------------------------------------------------------------
//                            fasmg_cycle
//------------------------------------------------------------------------

// This virtual table determines the implementation of the fasmg_cycle.
typedef struct
{
  // Executes a multigrid cycle using the various given algorithm components.
  void (*execute)(void* context, 
                  fasmg_operator_t* A, 
                  fasmg_grid_t* grid, 
                  fasmg_prolongator_t* prolongator,
                  fasmg_restrictor_t* restrictor,
                  real_t* B,
                  real_t* X);
  // Destructor.
  void (*dtor)(void* context);
} fasmg_cycle_vtable;

// Creates a fasmg_cycle that uses the given context and vtable to 
// implement its behavior. 
fasmg_cycle_t* fasmg_cycle_new(const char* name, 
                               void* context,
                               fasmg_cycle_vtable vtable);

// Frees the given fasmg_cycle.
void fasmg_cycle_free(fasmg_cycle_t* cycle);

// Returns the name of the cycle (internally stored).
char* fasmg_cycle_name(fasmg_cycle_t* cycle);

// Returns the context object for this integrator.
void* fasmg_cycle_context(fasmg_cycle_t* cycle);

// Performs a cycle on the given discretization for the given vectors using 
// the given operator/grid/prolongator/restrictor.
void fasmg_cycle_execute(fasmg_cycle_t* cycle,
                         fasmg_operator_t* A,
                         fasmg_grid_t* grid,
                         fasmg_prolongator_t* prolongator,
                         fasmg_restrictor_t* restrictor,
                         real_t* B,
                         real_t* X);

// Creates a new cycle object representing a V cycle with nu_1 relaxation steps 
// before the correction step, and nu_2 relaxation steps after.
fasmg_cycle_t* v_fasmg_cycle_new(int nu_1, int nu_2);

// Creates a new cycle object representing a W cycle with nu_1 relaxation steps
// before the correction step, and nu_2 relaxation steps after.
fasmg_cycle_t* w_fasmg_cycle_new(int nu_1, int nu_2);

// Creates a new cycle object representing a "mu" cycle with mu recursive 
// multigrid calls, nu_1 relaxation steps before the correction step, and 
// nu_2 relaxation steps after. In particular, mu = 1 is the V cycle and mu = 2 
// is the W cycle.
fasmg_cycle_t* mu_fasmg_cycle_new(int mu, int nu_1, int nu_2);

// Creates a new cycle object representing a full multi-grid (FMG) cycle that 
// executes nu_0 V cycles, each with nu_1 and nu_2 pre/post-correction 
// relaxation steps.
fasmg_cycle_t* fmg_fasmg_cycle_new(int nu_0, int nu_1, int nu_2);

#endif

