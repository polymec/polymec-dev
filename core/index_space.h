// Copyright (c) 2012-2013, Jeffrey N. Johnson
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

#ifndef POLYMEC_INDEX_SPACE_H
#define POLYMEC_INDEX_SPACE_H

#include "core/polymec.h"

// This data structure holds an index space that is distributed over a 
// set of parallel processes. Objects of this type are garbage-collected.
typedef struct 
{
  MPI_Comm comm; // Communicator on which index space is defined.
  int rank;      // Rank in the given communicator.
  int nproc;     // Number of parallel processes.
  int low, high; // Low and high indexes on this process.
} index_space_t;

// Creates an index space that has contiguous ranges on each parallel process
// for the given communicator, given the number of local indices. This 
// involves collective communication.
index_space_t* index_space_new(MPI_Comm comm, int num_local_indices);

// Creates an index space by naively partitioning the range [0, N] into 
// equal pieces, with the exception of the last process, which contains a 
// remainder.
index_space_t* index_space_from_naive_partitions(MPI_Comm comm, int N);

// Creates an index space from the given local data range.
index_space_t* index_space_from_low_and_high(MPI_Comm comm, int low, int high);

#endif

