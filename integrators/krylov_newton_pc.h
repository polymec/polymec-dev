// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_KRYLOV_NEWTON_PC_H
#define POLYMEC_KRYLOV_NEWTON_PC_H

#include "integrators/newton_pc.h"

// Types of Gram-Schmidt orthogonalization for GMRES preconditioners.
typedef enum
{
  GMRES_MODIFIED_GRAM_SCHMIDT,
  GMRES_CLASSICAL_GRAM_SCHMIDT
} gmres_newton_pc_gram_schmidt_t;

// The Krylov Newton PCs are truly matrix-free preconditioners. They can be 
// used with the Jacobian-Free Newton-Krylov, in particular the GMRES method, 
// which Knoll says is "surprisingly forgiving of mildly inconsistent 
// preconditioning." If you want to be more conservative, use Flexible GMRES
// (FGMRES) in your Newton-Krylov method.

// Creates a GMRES Newton preconditioner corresponding to the given function 
// F(t, x) OR the corresponding Jacobian J for which the product J*v is 
// specified. Exactly one of F and Jv must be non-NULL. If F is NULL, Jv is 
// used to compute the preconditioner matrix P. If Jv is NULL, a difference 
// quotient is used to compute J*v.
newton_pc_t* gmres_newton_pc_new(void* context,
                                 int (*F)(void* context, int i, real_t t, real_t* x, real_t* F),
                                 int (*Jv)(void* context, int i, real_t t, real_t* x, real_t* Jv),
                                 void (*dtor)(void* context),
                                 int num_local_values, 
                                 int num_remote_values,
                                 int max_restarts,
                                 gmres_newton_pc_gram_schmidt_t gs);
                                        
// Creates a version of the GMRES Newton preconditioner for solving DAEs 
// corresponding to the given function F(t, x, x_dot) = 0.
// OR the corresponding Jacobian J for which the product J*v is 
// specified. Exactly one of F and Jv must be non-NULL. If F is NULL, Jv is 
// used to compute the preconditioner matrix P. If Jv is NULL, a difference 
// quotient is used to compute J*v.
newton_pc_t* dae_gmres_newton_pc_new(void* context,
                                     int (*f)(void* context, int i, real_t t, real_t* x, real_t* x_dot, real_t* F),
                                     int (*Jv)(void* context, int i, real_t t, real_t* x, real_t* x_dot, real_t* Jv),
                                     void (*dtor)(void* context),
                                     int num_local_values, 
                                     int num_remote_values,
                                     int max_restarts,
                                     gmres_newton_pc_gram_schmidt_t gs);

bool newton_pc_is_gmres_newton_pc(newton_pc_t* gmres_pc);

void gmres_newton_pc_set_scaling(newton_pc_t* gmres_pc, real_t* scaling);
                                        
#endif

