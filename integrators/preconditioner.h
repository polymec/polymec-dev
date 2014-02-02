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

#ifndef POLYMEC_PRECONDITIONER_H
#define POLYMEC_PRECONDITIONER_H

#include "core/polymec.h"
#include "core/adj_graph.h"

// This class constructs and solves a preconditioned linear system to aid 
// the Jacobian-free Newton-Krylov solution of a nonlinear system.
typedef struct preconditioner_t preconditioner_t;

// This class holds information representing a preconditioner matrix.
typedef struct preconditioner_matrix_t preconditioner_matrix_t;

// This function constructs a new matrix for the given preconditioner.
typedef preconditioner_matrix_t* (*preconditioner_matrix_func)(void* context);

// This function computes the preconditioned Jacobian for the given preconditioner.
typedef void (*preconditioner_compute_jac_func)(void* context, real_t t, real_t* x, preconditioner_matrix_t* mat);

// This function solves the preconditioner system A*X = B, where A is the 
// preconditioner matrix and B is the given right-hand side. The right-hand 
// side is taken as input and is filled with the solution on output.
typedef void (*preconditioner_solve_func)(void* context, preconditioner_matrix_t* A, real_t* B);

// This function destroys the data context for the preconditioner.
typedef void (*preconditioner_dtor)(void* context);

// This virtual table determines the behavior of the preconditioner.
typedef struct
{
  preconditioner_matrix_func matrix;
  preconditioner_compute_jac_func compute_jacobian;
  preconditioner_solve_func solve;
  preconditioner_dtor dtor;
} preconditioner_vtable;

// Constructs a preconditioner with the given name, data context, 
// virtual table. 
preconditioner_t* preconditioner_new(const char* name,
                                     void* context,
                                     preconditioner_vtable vtable);

// Destroys the given preconditioner.
void preconditioner_free(preconditioner_t* precond);

// Returns an internal pointer to the name of the preconditioner.
char* preconditioner_name(preconditioner_t* precond);

// Returns the data context for the preconditioner.
void* preconditioner_context(preconditioner_t* precond);

// Constructs a new preconditioner matrix for this preconditioner.
preconditioner_matrix_t* preconditioner_matrix(preconditioner_t* precond);

// Calculates the coefficients of the preconditioner matrix corresponding to 
// the Jacobian of the underlying nonlinear system, storing the result in mat.
void preconditioner_compute_jacobian(preconditioner_t* precond,
                                     real_t t,
                                     real_t* x,
                                     preconditioner_matrix_t* mat);

// Solves the preconditioned linear system. On input, rhs contains the 
// right-hand side of the system, and on output, it contains the solution 
// to the preconditioned system.
void preconditioner_solve(preconditioner_t* precond,
                          preconditioner_matrix_t* mat,
                          real_t* rhs);

// Transforms the given preconditioner matrix A to (I - gamma * A).
typedef void (*preconditioner_matrix_scale_and_shift_func)(void* context, real_t gamma);

// This function provides (read-only) access to the (i, j)th coefficient 
// in a preconditioner matrix.
typedef real_t (*preconditioner_matrix_coeff_func)(void* context, int i, int j);

// This function destroys the data context for the preconditioner matrix.
typedef void (*preconditioner_matrix_dtor)(void* context);

// This virtual table determines the nature of the preconditioner matrix 
// representation.
typedef struct
{
  preconditioner_matrix_scale_and_shift_func scale_and_shift;
  preconditioner_matrix_coeff_func           coeff;
  preconditioner_matrix_dtor                 dtor;
} preconditioner_matrix_vtable;

// Constructs a representation of a (square) preconditioner matrix with the 
// given name, data context, virtual table, and number of rows. The representation
// of the matrix is tied to the implementation of a particular preconditioner.
preconditioner_matrix_t* preconditioner_matrix_new(const char* name,
                                                   void* context,
                                                   preconditioner_matrix_vtable vtable,
                                                   int num_rows);

// Destroys the given preconditioner matrix.
void preconditioner_matrix_free(preconditioner_matrix_t* mat);

// Returns the data context for the preconditioner matrix.
void* preconditioner_matrix_context(preconditioner_matrix_t* mat);

// Transforms the given preconditioner matrix A to (I - gamma * A), where 
// I is the identity matrix. This is used for preconditioning time-dependent
// problems.
void preconditioner_matrix_scale_and_shift(preconditioner_matrix_t* mat, real_t gamma);

// Returns the (i, j)th entry in the preconditioner matrix. This method of 
// access is slow in general and should only be used for diagnostics.
real_t preconditioner_matrix_coeff(preconditioner_matrix_t* mat, int i, int j);

#endif

