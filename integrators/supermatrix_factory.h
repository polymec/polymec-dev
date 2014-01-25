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

#ifndef POLYMEC_SUPERMATRIX_FACTORY_H
#define POLYMEC_SUPERMATRIX_FACTORY_H

#include "core/polymec.h"
#include "core/adj_graph.h"
#include "cvode/cvode.h"
#include "slu_ddefs.h"
#include "supermatrix.h"

// This class represents a way of integrating a system of 
// nonlinear equations.
typedef struct supermatrix_factory_t supermatrix_factory_t;

// Creates a factory that produces SuperLU SuperMatrix objects using 
// the given adjacency graph and system function F (and context). The 
// factory does NOT assume ownership of the factory, graph, or context.
supermatrix_factory_t* supermatrix_factory_new(adj_graph_t* graph,
                                               int (*sys_func)(void* context, real_t t, real_t* x, real_t* F), 
                                               void* context);

// Frees the factory.
void supermatrix_factory_free(supermatrix_factory_t* factory);

// Produces a supermatrix with a sparsity that matches the given graph 
// from which the factory was constructed.
SuperMatrix* supermatrix_factory_matrix(supermatrix_factory_t* factory);

SuperMatrix* supermatrix_factory_vector(supermatrix_factory_t* factory,
                                        int num_rhs);

// Produces a Jacobian supermatrix from finite difference quotients applied 
// to the given function F at the solution x and time t.
SuperMatrix* supermatrix_factory_jacobian(supermatrix_factory_t* factory, real_t* x, real_t t);

// Updates a Jacobian supermatrix from finite difference quotients applied 
// to the given function F at the solution x and time t.
void supermatrix_factory_update_jacobian(supermatrix_factory_t* factory, real_t* x, real_t t, SuperMatrix* J);

// Call this function to destroy supermatrices that have been created by this 
// factory.
void supermatrix_free(SuperMatrix* matrix);

#endif

