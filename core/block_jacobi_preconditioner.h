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

#ifndef POLYMEC_BLOCK_JACOBI_PRECONDITIONER_H
#define POLYMEC_BLOCK_JACOBI_PRECONDITIONER_H

#include "core/preconditioner.h"
#include "core/adj_graph.h"

// Creates a block Jacobi preconditioner from a rule that computes the 
// (block) diagonal. The diagonal D is stored in block-major order. Specifically, 
// it is an array consisting of the components of column-major dense matrices 
// stacked end-to-end.
preconditioner_t* block_jacobi_preconditioner_new(void* context,
                                                  void (*compute_diagonal)(void* context, int block_size, real_t alpha, real_t beta, real_t gamma, real_t t, real_t* x, real_t* x_dot, real_t* D),
                                                  void (*dtor)(void* context),
                                                  int num_block_rows,
                                                  int block_size);

// Use this to retrieve the context pointer for the given block Jacobi 
// preconditioner. DON'T use this unless you know what you're doing.
void* block_jacobi_preconditioner_context(preconditioner_t* bj_precond);

// Creates a block Jacobi preconditioner with the given sparsity graph, 
// number of block rows, and block size. The nature of the sparsity graph 
// (i.e. whether it is a block graph or not) is inferred from the number of 
// block rows and the block size.
preconditioner_t* block_jacobi_preconditioner_from_function(void* context,
                                                            int (*F)(void* context, real_t t, real_t* x, real_t* Fval),
                                                            void (*dtor)(void* context),
                                                            adj_graph_t* sparsity,
                                                            int num_block_rows,
                                                            int block_size);
                                        
// Creates a block Jacobi preconditioner with the given sparsity graph, 
// number of block rows, and block size, appropriate for preconditioning 
// differential-algebraic systems.
preconditioner_t* block_jacobi_dae_preconditioner_from_function(void* context,
                                                                int (*F)(void* context, real_t t, real_t* x, real_t* xdot, real_t* Fval),
                                                                void (*dtor)(void* context),
                                                                adj_graph_t* sparsity,
                                                                int num_block_rows,
                                                                int block_size);
                                        
#endif

