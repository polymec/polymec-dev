#ifndef POLYMEC_INDEX_SPACE_H
#define POLYMEC_INDEX_SPACE_H

#include "core/polymec.h"

#ifdef __cplusplus
extern "C" {
#endif

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

#ifdef __cplusplus
}
#endif

#endif

