// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_SSOR_NEWTON_PC_H
#define POLYMEC_SSOR_NEWTON_PC_H

#include "integrators/newton_pc.h"

// The SSOR Newton PC implements a nonlinear Symmetric Successive 
// Overrelaxation Newton preconditioning method described by Chan and 
// Jackson, Siam J. Sci. Stat. Comput. 5 (1984). It is mostly-matrix-free
// preconditioner in that it requires only the diagonal of the Jacobian to 
// be computed.

// Creates a SSOR Newton preconditioner with the given relaxation factor 
// omega, f computes the value of the nonlinear function f(t, x) = 0 at a given 
// index i, DJ computes the diagonal elements of the Jacobian of f at (t, x)
// at index i, and converged returns true if the iteration has converged at 
// index i and false if not.
newton_pc_t* ssor_newton_pc_new(void* context,
                                real_t (*f)(void* context, int i, real_t t, real_t* x),
                                real_t (*DJ)(void* context, int i, real_t t, real_t* x),
                                void (*dtor)(void* context),
                                int num_local_values,
                                real_t omega);
                                        
// This function returns true if the given Newton preconditioner object 
// is an SSOR Newton preconditioner, false if not.
bool newton_pc_is_ssor_newton_pc(newton_pc_t* newton_pc);

// This function sets the relaxation factor to the given value omega within 
// the range (0, 2), provided that the given Newton preconditioner is an SSOR 
// preconditioner.
void ssor_newton_pc_set_relaxation_factor(newton_pc_t* ssor_newton_pc,
                                          real_t omega);

#endif

