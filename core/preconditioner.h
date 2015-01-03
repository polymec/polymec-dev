// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_PRECONDITIONER_H
#define POLYMEC_PRECONDITIONER_H

#include "core/polymec.h"

// This class constructs and solves a preconditioned linear system to aid 
// the Jacobian-free Newton-Krylov solution of a nonlinear system.
typedef struct preconditioner_t preconditioner_t;

// This function sets up the preconditioner to obtain a preconditioned solution.
typedef void (*preconditioner_setup_func)(void* context);

// This function solves the preconditioner system A*X = B, where A is the 
// preconditioner matrix and B is the given right-hand side. The right-hand 
// side is taken as input and is filled with the solution on output.
typedef bool (*preconditioner_solve_func)(void* context, real_t* B);

// This function writes a text representation of the preconditioner to the given file.
typedef void (*preconditioner_fprintf_func)(void* context, FILE* stream);

// This function destroys the data context for the preconditioner.
typedef void (*preconditioner_dtor)(void* context);

// This virtual table determines the behavior of the preconditioner.
typedef struct
{
  preconditioner_setup_func   setup;
  preconditioner_solve_func   solve;
  preconditioner_fprintf_func fprintf;
  preconditioner_dtor         dtor;
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

// Sets up the preconditioner to obtain a representation of the preconditioner
// operator. 
void preconditioner_setup(preconditioner_t* precond);

// Solves the preconditioned linear system. On input, rhs contains the 
// right-hand side of the system, and on output, it contains the solution 
// to the preconditioned system. Returns true if the system was successfully 
// solved, false if not.
bool preconditioner_solve(preconditioner_t* precond, real_t* rhs);

// Writes a text representation of the given preconditioner to the given file.
void preconditioner_fprintf(preconditioner_t* precond, FILE* stream);

#endif

