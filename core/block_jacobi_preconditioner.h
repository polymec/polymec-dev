// Copyright (c) 2012-2015, Jeffrey N. Johnson
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
                                                  void (*compute_diagonal)(void* context, int block_size, real_t* D),
                                                  void (*dtor)(void* context),
                                                  int num_block_rows,
                                                  int block_size);

// Creates a block Jacobi preconditioner that uses data in the given 
// (block-major) array as its diagonal. The preconditioner does NOT manage the 
// array -- this array is assumed to be valid for the lifetime of the preconditioner.
preconditioner_t* block_jacobi_preconditioner_from_array(real_t* array,
                                                         int num_block_rows,
                                                         int block_size);

#endif

