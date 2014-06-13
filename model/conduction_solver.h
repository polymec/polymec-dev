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

#ifndef POLYMEC_CONDUCTION_SOLVER_H
#define POLYMEC_CONDUCTION_SOLVER_H

#include "core/krylov_solver.h"
#include "core/mesh.h"

// The conduction solver is a (Krylov) linear solver that uses a low-order 
// finite-volume discretization to obtain a solution to the steady-state 
// heat equation 
// 
//                          div (lambda o grad phi) = f 
//
// on a mesh. The solution phi and the quantities lambda and f are assumed to be 
// cell-centered scalars. This solver is a fine way to produce good initial 
// guesses for nonlinear elliptical partial differential equations that 
// resemble the steady-state heat equation (or any of its variants).

// These constructors all create a new conduction solver that solve the 
// steady-state heat equation on the given mesh. 

// Creates a solver that uses the GMRES method with the given max Krylov 
// subspace dimension and the maximum number of restarts.
krylov_solver_t* gmres_conduction_solver_new(mesh_t* mesh, int max_krylov_dim, int max_restarts);

// Creates a solver that uses the BiCGSTAB method with the given max Krylov 
// subspace dimension.
krylov_solver_t* bicgstab_conduction_solver_new(mesh_t* mesh, int max_krylov_dim);

// Creates a solver that uses the TFQMR method with the given max Krylov 
// subspace dimension.
krylov_solver_t* tfqmr_conduction_solver_new(mesh_t* mesh, int max_krylov_dim);

// Returns an internal pointer to the array that stores the cell-centered 
// values of lambda. By default, lambda is set to 1 everywhere.
real_t* conduction_solver_lambda(krylov_solver_t* solver);

// Returns an internal pointer to the array that stores the cell-centered 
// values of f. By default, f is set to 0 everywhere.
real_t* conduction_solver_f(krylov_solver_t* solver);

// Adds a Robin-type boundary condition to the faces in the solver's mesh 
// corresponding to the given tag: 
//   alpha * phi + beta * n o grad phi = gamma
void conduction_solver_add_bc(krylov_solver_t* solver, const char* face_tag);

// Returns internal pointers to the alpha/beta/gamma arrays that define the 
// values of these coefficients for the boundary condition assigned to the 
// given face tag. Also stores the number of faces in the tag in num_faces.
// The arrays alpha, beta, and gamma store values for the faces in the given 
// face tag, such that (e. g.) alpha[i] stores the value of alpha on the face 
// with index tag[i] if tag is the array associated with face_tag. 
// By default, the values in all arrays are set to zero. Any of these array 
// pointers can be set to NULL if no data is requested. If no such tag is 
// found, these array pointer values are set to NULL.
void conduction_solver_get_bc_arrays(krylov_solver_t* solver, const char* face_tag, 
                                     real_t** alpha, real_t** beta, real_t** gamma, int* num_faces);

#endif
