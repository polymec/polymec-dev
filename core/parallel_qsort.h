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

#ifndef POLYMEC_PARALLEL_QSORT_H
#define POLYMEC_PARALLEL_QSORT_H

#include "core/rng.h"

// This function implements a parallel quicksort using regular or 
// random sampling. Data is assumed to exist for each process belonging to 
// the communicator comm, in an array pointed to by base, each having a 
// (locally-stored) number of elements nel, each of a size width in bytes. 
// The comparator compar is used to perform comparisons between local elements 
// in each array. If rng is NULL, sampling for pivots is performed at regular 
// intervals. Otherwise, rng is used to randomly select pivots on each process.
// On output, base points to the same locally-stored array as before, but
// its data has been sorted on each process p such that processes preceding p 
// contain sorted data preceding that on p, and processes following p contain 
// sorted data following that on p.
void parallel_qsort(MPI_Comm comm, 
                    void* base, 
                    size_t nel, 
                    size_t width,
                    int (*compar)(const void* left, const void* right),
                    rng_t* rng);

#endif

