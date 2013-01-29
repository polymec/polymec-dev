#include <gc/gc.h>
#include "core/index_space.h"

#ifdef __cplusplus
extern "C" {
#endif

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

#ifdef __cplusplus
}
#endif

