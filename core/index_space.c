#include "core/index_space.h"

#ifdef __cplusplus
extern "C" {
#endif

index_space_t* index_space_new(MPI_Comm comm, int num_local_indices)
{
  int rank, nproc;
  MPI_Comm_size(comm, &nproc);
  MPI_Comm_rank(comm, &rank);
  int sizes[nproc];
  MPI_Allgather(&num_local_indices, 1, MPI_INT, sizes, 1, MPI_INT, comm);

  index_space_t* space = malloc(sizeof(index_space_t));
  space->comm = comm;
  space->low = 0;
  for (int i = 0; i < rank; ++i)
    space->low += sizes[i];
  space->high = space->low + num_local_indices;
  return space;
}

void index_space_free(index_space_t* index_space)
{
  free(index_space);
}

#ifdef __cplusplus
}
#endif

