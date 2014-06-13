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
// steady-state heat equation on the given mesh. The arrays lambda and f are 
// used by the solver to provide the values of lambda and f on cell centers, 
// but are not consumed by the solver. If lambda is NULL, Poisson's equation 
// is solved. If f is NULL, the homogeneous steady-state heat equation is 
// solved.

// Creates a solver that uses the GMRES method with the given max Krylov 
// subspace dimension and the maximum number of restarts.
krylov_solver_t* gmres_conduction_solver_new(mesh_t* mesh, real_t* lambda, real_t* f,
                                             int max_krylov_dim, int max_restarts);

// Creates a solver that uses the BiCGSTAB method with the given max Krylov 
// subspace dimension.
krylov_solver_t* bicgstab_conduction_solver_new(mesh_t* mesh, real_t* lambda, real_t* f,
                                                int max_krylov_dim);

// Creates a solver that uses the TFQMR method with the given max Krylov 
// subspace dimension.
krylov_solver_t* tfqmr_conduction_solver_new(mesh_t* mesh, real_t* lambda, real_t* f,
                                            int max_krylov_dim);

// Adds a Robin-type boundary condition to the faces in the solver's mesh 
// corresponding to the given tag: 
//
// alpha * phi + beta * n o grad phi = gamma
//
// The arrays alpha, beta, and gamma store values for the faces in the given 
// face tag, such that (e. g.) alpha[i] stores the value of alpha on the face 
// with index tag[i] if tag is the array associated with face_tag. The arrays 
// alpha, beta, and gamma are used by the solver but not consumed by it.
void conduction_solver_add_bc(krylov_solver_t* conduction, const char* face_tag,
                              real_t* alpha, real_t* beta, real_t* gamma);


#endif
