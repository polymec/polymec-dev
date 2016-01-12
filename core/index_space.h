// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

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

