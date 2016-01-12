// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_PHYSICS_KERNEL_H
#define POLYMEC_PHYSICS_KERNEL_H

#include "core/krylov_solver.h"
#include "multiphysics/physics_state.h"

// This class provides an abstract interface for a physics kernel that is 
// responsible for evolving the state of a system according to the physics
// for a single process.
typedef struct physics_kernel_t physics_kernel_t;

// A physics kernel can be one of three types, depending on the nature of 
// the underlying process, and the corresponding time scale: 
// (1) A HYPERBOLIC kernel represents a "slow" process that occurs at the 
//     time scale of the integration. Hyperbolic processes are so named for 
//     their complex (plane-wave) eigenvalues and represent wave-like motions.
// (2) A PARABOLIC kernel represents a faster process or set of processes that 
//     occur at one or more frequencies, ranging from the time scale of 
//     integration to much faster, but still clearly evolving continously. 
//     Parabolic processes look like diffusion.
// (3) An ELLIPTIC kernel represents a very fast process that achieves a 
//     steady state instantaneously at the time scale of integration. 
//     Elliptic processes behave like relaxation process and/or constraints.
typedef enum
{
  PHYSICS_KERNEL_HYPERBOLIC,
  PHYSICS_KERNEL_PARABOLIC,
  PHYSICS_KERNEL_ELLIPTIC
} physics_kernel_type_t;

// This virtual table determines the implementation of the kernel.
typedef struct
{
  // This function returns the size of the given primary variable for this kernel.
  int (*primary_size)(void* context, const char* var_name);

  // This function returns the size of the given secondary variable for this kernel.
  int (*secondary_size)(void* context, const char* var_name);

  // This function computes the time change in the given array of primary 
  // variables u and the secondary variables w, ordered as defined during its setup, 
  // at the time t. This must be provided for HYPERBOLIC or PARABOLIC physics kernels.
  void (*eval_dudt)(void* context, real_t t, real_t* u, real_t* w, real_t* dudt);

  // This function computes the residual function R = F(t, u) ~ 0 that 
  // characterizes the constraint/relaxation process represented by the 
  // physics kernel. Primary variables are given in the vector u, while secondaries are
  // given in w. It must be provided for ELLIPTIC physics kernels.
  void (*eval_residual)(void* context, real_t t, real_t* u, real_t* w, real_t* R);

  // This function returns the maximum timestep that can be used to integrate 
  // the physical process represented by this kernel, given primary variables u and 
  // secondary variables w at time t. Must be provided for HYPERBOLIC physics kernels.
  real_t (*max_dt)(void* context, real_t t, real_t* u, real_t* w);

  //------------------------------------------------------------------------
  //                      Jacobian-related functions
  //------------------------------------------------------------------------
  // What is meant by the Jacobian depends on whether the physics kernel is 
  // HYPERBOLIC/PARABOLIC or ELLIPTIC. For HYPERBOLIC/PARABOLIC kernels, 
  // the Jacobian is the derivative of the time derivative du/dt w.r.t. u.
  // For ELLIPTIC kernels, the Jacobian is the derivative of the residual 
  // function R = F(t, u) ~ 0 w.r.t. u.
  //------------------------------------------------------------------------

  // This optional function computes the matrix-vector product of the 
  // Jacobian with a vector v, given the primary variables u and secondary variables 
  // w at time t, placing the result in Jv. If it is not provided, a differencing 
  // approximation will be used.
  void (*compute_Jv)(void* context, real_t t, real_t* u, real_t* w, real_t* v, real_t* Jv);

  // This optional function computes the components of the Jacobian 
  // at time t for the primary variables u and secondary variables w, placing them into 
  // the given Krylov matrix. The matrix must have been created with the proper sparsity 
  // pattern.
  void (*compute_J)(void* context, real_t t, real_t* u, real_t* w, krylov_matrix_t* J);

  // This function destroys the state (context) when the kernel 
  // is destroyed.
  void (*dtor)(void* context);

} physics_kernel_vtable;

// This type represents a function that is used with a physics kernel's context to 
// update a given secondary variable at time t, given primary and secondary variable 
// vectors u and w.
typedef void (*physics_kernel_update_function_t)(void* context, real_t t, real_t* u, real_t* w, real_t* secondary_var);

// Creates a physics kernel of the given type that uses the given context and 
// vtable.
physics_kernel_t* physics_kernel_new(const char* name, 
                                     physics_kernel_type_t type,
                                     void* context,
                                     physics_kernel_vtable vtable);

// Frees a physics kernel.
void physics_kernel_free(physics_kernel_t* kernel);

//------------------------------------------------------------------------
//               Basic attributes of the physics kernel
//------------------------------------------------------------------------

// Returns the name of the physics kernel (internally stored).
char* physics_kernel_name(physics_kernel_t* kernel);

// Returns the type of the physics kernel.
physics_kernel_type_t physics_kernel_type(physics_kernel_t* kernel);

// Returns the context object for this kernel.
void* physics_kernel_context(physics_kernel_t* kernel);

// Returns the size of the solution vector for this kernel.
int physics_kernel_solution_vector_size(physics_kernel_t* kernel);

//------------------------------------------------------------------------
//                  Construction metadata interface
//------------------------------------------------------------------------
// The following methods should be called only during the construction of 
// a physics kernel.
//------------------------------------------------------------------------

// Adds a primary variable with the given number of components to the list 
// of evolved primary variables for the physics kernel. Indexing within 
// the primary variable vector starts from 0 and is incremented for each 
// primary variable added.
void physics_kernel_add_primary(physics_kernel_t* kernel,
                                const char* var_name,
                                int num_components);

// Returns true if the physics kernel uses a primary variable with the given name, 
// false otherwise.
bool physics_kernel_has_primary(physics_kernel_t* kernel, const char* var_name);

// Returns the size of the primary variable with the given name within the 
// kernel, or -1 if the kernel contains no such variable.
int physics_kernel_primary_size(physics_kernel_t* kernel, const char* var_name);

// Traverses the list of this physics kernel's primary variables, retrieving 
// each one's name, index within the primary vector u, its size (number of data)
// and its number of components per datum. Returns true if the traversal 
// yielded another primary variable, false if not. *pos must be set to zero 
// for the traversal to be reset.
bool physics_kernel_next_primary(physics_kernel_t* kernel,
                                 int* pos,
                                 char** var_name,
                                 int* index,
                                 int* size,
                                 int* num_components);

// Adds a secondary variable with the given number of components to the list 
// of maintained secondary variables for the physics kernel. Uses the given 
// update function to update the secondary variable as needed. This update 
// function uses the context for the kernel to perform its work, as well 
// as the primary vector u and the time t.
void physics_kernel_add_secondary(physics_kernel_t* kernel,
                                  const char* var_name,
                                  int num_components,
                                  physics_kernel_update_function_t update);

// Adds a (secondary variable) dependency to the given secondary variable for 
// the physics kernel. That is: in order for the secondary variable var_name 
// to be updated, the secondary variable dep_name must be updated first. It is
// an error for either of these secondary variables not to exist in the kernel.
void physics_kernel_add_secondary_dep(physics_kernel_t* kernel,
                                      const char* var_name,
                                      const char* dep_name);

// Returns true if the physics kernel uses and/or maintains a secondary variable with 
// the given name, false otherwise.
bool physics_kernel_has_secondary(physics_kernel_t* kernel, const char* var_name);

// Returns the size of the secondary variable with the given name within the 
// kernel, or -1 if the kernel contains no such variable.
int physics_kernel_secondary_size(physics_kernel_t* kernel, const char* var_name);

// Traverses the list of this physics kernel's secondary variables, retrieving 
// each one's name, index, number of data, number of components per data, and 
// update function. Returns true if the traversal yielded another secondary 
// variable, false if not. *pos must be set to zero for the traversal to be 
// reset.
bool physics_kernel_next_secondary(physics_kernel_t* kernel,
                                   int* pos,
                                   char** var_name,
                                   int* index,
                                   int* size,
                                   int* num_components,
                                   physics_kernel_update_function_t* update);

// Traverses the list of dependencies for the given secondary variable in 
// this physics kernel, returning true if a dependency was yielded, false if 
// not. *pos must be set to zero for the traversal to be reset.
bool physics_kernel_next_secondary_dep(physics_kernel_t* kernel,
                                       const char* var_name,
                                       int* pos,
                                       char** dep_name);

// Updates the given secondary variable data within the given state using 
// the update function registered within this kernel for this variable, 
// and at the time t. It is an error for this variable not to exist in this 
// kernel.
void physics_kernel_update_secondary(physics_kernel_t* kernel,
                                     const char* var_name,
                                     real_t t,
                                     physics_state_t* state);

//------------------------------------------------------------------------
//                      Time integration interface
//------------------------------------------------------------------------
// The following methods are called to interact with a multiphysics 
// time integrator.
//------------------------------------------------------------------------

// Returns the current number of data for the given primary variable.
int physics_kernel_primary_size(physics_kernel_t* kernel, const char* var_name);

// Returns the current number of data for the given secondary variable.
int physics_kernel_secondary_size(physics_kernel_t* kernel, const char* var_name);

// Returns the block size of the Jacobian matrix associated with this physics
// kernel, which should be the sum of all of the components of the primary 
// variables.
int physics_kernel_jacobian_block_size(physics_kernel_t* kernel);

// Computes the time derivative of the primary variable vector u, given the state at 
// time t. This may only be called on HYPERBOLIC or PARABOLIC physics kernels.
void physics_kernel_eval_dudt(physics_kernel_t* kernel, real_t t, physics_state_t* state, real_t* dudt);

// Computes the value of the residual function R = F(t, u) for the state at time t. This 
// may only be called on ELLIPTIC physics kernels.
void physics_kernel_eval_residual(physics_kernel_t* kernel, real_t t, physics_state_t* state, real_t* R);

// Returns the maximum time step size that can be used to advance the 
// solution using this kernel for the given state at the given time t, 
// or FLT_MAX if no such maximum exists.
real_t physics_kernel_max_dt(physics_kernel_t* kernel, real_t t, physics_state_t* state);

// Computes the matrix-vector product of the Jacobian with a vector v at 
// the given time t for the given state.
void physics_kernel_compute_Jv(physics_kernel_t* kernel, real_t t, physics_state_t* state, real_t* v, real_t* Jv);

// Computes the components of the Jacobian at time t for the given state,
// placing them into the given Krylov matrix. The matrix must have been 
// created with the proper sparsity pattern.
void physics_kernel_compute_J(physics_kernel_t* kernel, real_t t, physics_state_t* state, krylov_matrix_t* J);

#endif

