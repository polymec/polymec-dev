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

#include <gc/gc.h>
#include "core/index_space.h"

index_space_t* index_space_new(MPI_Comm comm, int num_local_indices)
{
#if USE_MPI
  int rank, nproc;
  MPI_Comm_size(comm, &nproc);
  MPI_Comm_rank(comm, &rank);
  int sizes[nproc];
  MPI_Allgather(&num_local_indices, 1, MPI_INT, sizes, 1, MPI_INT, comm);
  int low = 0, high;
  for (int i = 0; i < rank; ++i)
    low += sizes[i];
  high = space->low + num_local_indices;
#else
  int nproc = 1, rank = 0, low = 0, high = num_local_indices;
#endif

  index_space_t* space = GC_MALLOC(sizeof(index_space_t));
  space->comm = comm;
  space->rank = rank;
  space->nproc = nproc;
  space->low = low;
  space->high = high;
  return space;
}

index_space_t* index_space_from_naive_partitions(MPI_Comm comm, int N)
{
#if USE_MPI
  int rank, nproc;
  MPI_Comm_size(comm, &nproc);
  MPI_Comm_rank(comm, &rank);
  int ipp = N / nproc;
  int low = rank * ipp;
  int high = low + ipp; 
  if (rank == (nproc - 1))
    high = N;
#else
  int nproc = 1, rank = 0, low = 0, high = N;
#endif

  index_space_t* space = GC_MALLOC(sizeof(index_space_t));
  space->comm = comm;
  space->rank = rank;
  space->nproc = nproc;
  space->low = low;
  space->high = high;
  return space;
}

index_space_t* index_space_from_low_and_high(MPI_Comm comm, int low, int high)
{
#if USE_MPI
  int rank, nproc;
  MPI_Comm_size(comm, &nproc);
  MPI_Comm_rank(comm, &rank);
#else
  int nproc = 1, rank = 0;
#endif

  index_space_t* space = GC_MALLOC(sizeof(index_space_t));
  space->comm = comm;
  space->rank = rank;
  space->nproc = nproc;
  space->low = low;
  space->high = high;
  return space;
}

