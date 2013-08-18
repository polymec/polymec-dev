// Copyright 2012-2013 Jeffrey Johnson.
// 
// This file is part of Polymec, and is licensed under the Apache License, 
// Version 2.0 (the "License"); you may not use this file except in 
// compliance with the License. You may may find the text of the license in 
// the LICENSE file at the top-level source directory, or obtain a copy of 
// it at
// 
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

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

