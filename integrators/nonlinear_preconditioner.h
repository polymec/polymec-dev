// Copyright (c) 2012-2014, Jeffrey N. Johnson
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this 
// list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice, 
// this list of conditions and the following disclaimer in the documentation 
// and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef POLYMEC_NONLINEAR_PRECONDITIONER_H
#define POLYMEC_NONLINEAR_PRECONDITIONER_H

#include "core/polymec.h"
#include "core/adj_graph.h"
#include "core/preconditioner.h"

// Every nonlinear preconditioner has to implement the functions in the 
// following virtual table.
typedef struct
{
  // Method for computing the preconditioner operator 
  // P = alpha * I + beta * dF/dx + gamma * dF/dxdot.
  void (*compute_P)(void* context, real_t alpha, real_t beta, real_t gamma, real_t t, real_t* x, real_t* xdot);
  // Method for solving the preconditioner system P * z = r. On input, z is 
  // the right-hand side, and on output it is the solution.
  bool (*solve)(void* context, real_t* z);
  // fprintf function.
  void (*fprintf)(void* context, FILE* stream);
  // Context destructor.
  void (*dtor)(void* context);
} nonlinear_preconditioner_vtable;

// Constructs a preconditioner that represents the linear combination
// alpha * I + beta * dF/dx + gamma * dF/dx, where F(t, x[, xdot]) = 0 is a 
// function representing a nonlinear ODE [or differential/algebraic) system, 
// I is the identity matrix, and dF/dx [and dF/dxdot] are the derivatives of 
// F with respect to x [and xdot].
preconditioner_t* nonlinear_preconditioner_new(const char* name,
                                               void* context,
                                               nonlinear_preconditioner_vtable vtable);

// Sets up the preconditioner at the point (t, x, xdot) in solution space, computing
// alpha * I + beta * dF/dx + gamma * dF/dx for the given alpha, beta, and gamma.
// Note that xdot can be NULL if F is only a function of x.
void nonlinear_preconditioner_setup(preconditioner_t* precond, 
                                    real_t alpha, real_t beta, real_t gamma,
                                    real_t t, real_t* x, real_t* xdot);

//------------------------------------------------------------------------
//                  Nonlinear preconditioner implementations
//------------------------------------------------------------------------

// Creates a block Jacobi preconditioner with the given sparsity graph, 
// number of block rows, and block size. The nature of the sparsity graph 
// (i.e. whether it is a block graph or not) is inferred from the number of 
// block rows and the block size.
preconditioner_t* block_jacobi_preconditioner_from_function(const char* name, 
                                                            void* context,
                                                            int (*F)(void* context, real_t t, real_t* x, real_t* Fval),
                                                            void (*dtor)(void* context),
                                                            adj_graph_t* sparsity,
                                                            int num_block_rows,
                                                            int block_size);
                                        
// Creates a block Jacobi preconditioner with the given sparsity graph, 
// number of block rows, and block size, appropriate for preconditioning 
// differential-algebraic systems.
preconditioner_t* block_jacobi_preconditioner_from_dae_function(const char* name, 
                                                                void* context,
                                                                int (*F)(void* context, real_t t, real_t* x, real_t* xdot, real_t* Fval),
                                                                void (*dtor)(void* context),
                                                                adj_graph_t* sparsity,
                                                                int num_block_rows,
                                                                int block_size);

// Sparse (Supernode) LU preconditioner.
preconditioner_t* lu_preconditioner_from_function(const char* name, 
                                                  void* context,
                                                  int (*F)(void* context, real_t t, real_t* x, real_t* Fval),
                                                  void (*dtor)(void* context),
                                                  adj_graph_t* sparsity,
                                                  int num_block_rows,
                                                  int block_size);
 
// Sparse (Supernode) LU preconditioner -- differential-algebraic version.
preconditioner_t* lu_preconditioner_from_dae_function(const char* name,
                                                      void* context,
                                                      int (*F)(void* context, real_t t, real_t* x, real_t* xdot, real_t* Fval),
                                                      void (*dtor)(void* context),
                                                      adj_graph_t* sparsity,
                                                      int num_block_rows,
                                                      int block_size);
 
// The following types give options to control ILU preconditioners for the 
// nonlinear and time integrators.

// Specifies how Incomplete LU (ILU) solvers should permute the rows of the 
// matrix.
typedef enum
{
  ILU_NO_ROW_PERM,      // no row permutations
  ILU_LARGE_DIAG_PERM   // try to make the diagonal larger
} ilu_row_perm_t;

// Species the type of (L-) norm used to measure errors in ILU solves.
typedef enum
{
  ILU_L1,  
  ILU_L2,  
  ILU_LINF
} ilu_norm_t;

// Rules for dropping coefficients in ILU. See SuperLU documentation.
// These rules can be combined with bitwise OR.
extern const int ILU_DROP_BASIC;
extern const int ILU_DROP_PROWS;
extern const int ILU_DROP_COLUMN;
extern const int ILU_DROP_AREA;
extern const int ILU_DROP_DYNAMIC;
extern const int ILU_DROP_INTERP;

// Identifies a variant (if any) of modified ILU.
typedef enum
{
  ILU_SILU,           // no modified ILU
  ILU_MILU1,
  ILU_MILU2,
  ILU_MILU3
} ilu_variant_t;

// This container specifies a complete set of ILU parameters. Objects 
// of this type are garbage-collected.
typedef struct
{
  real_t diag_pivot_threshold;
  ilu_row_perm_t row_perm;
  int drop_rule;
  real_t drop_tolerance;
  real_t fill_factor;
  ilu_variant_t milu_variant;
  real_t fill_tolerance;
  ilu_norm_t norm;
} ilu_params_t;

// Initializes a set of ILU parameters with reasonable defaults.
ilu_params_t* ilu_params_new();

// ILU preconditioner.
preconditioner_t* ilu_preconditioner_from_function(const char* name,
                                                   void* context,
                                                   int (*F)(void* context, real_t t, real_t* x, real_t* Fval),
                                                   void (*dtor)(void* context),
                                                   adj_graph_t* sparsity, 
                                                   int num_block_rows,
                                                   int block_size,
                                                   ilu_params_t* ilu_params);

// ILU preconditioner -- differential-algebraic version.
preconditioner_t* ilu_preconditioner_from_dae_function(const char* name,
                                                       void* context,
                                                       int (*F)(void* context, real_t t, real_t* x, real_t* xdot, real_t* Fval),
                                                       void (*dtor)(void* context),
                                                       adj_graph_t* sparsity, 
                                                       int num_block_rows,
                                                       int block_size,
                                                       ilu_params_t* ilu_params);

#endif

