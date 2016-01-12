// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

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

