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

#include "core/parallel_qsort_with_sampling.h"

// This helper merges an array of sorted lists of data (having elements of 
// the given width in bytes), merging them into a single sorted list using 
// the given comparator function.
static void merge_sorted_lists(uint8_t** lists_to_merge, 
                               int* list_sizes,
                               int num_lists,
                               int width,
                               int (*compar)(const void* left, const void* right),
                               uint8_t* merged_list)
{
  // For now, we do the incredibly stupid thing--just make a giant list 
  // out of all the smaller ones, and then qsort it.
  // FIXME: Suboptimal much??
  int k = 0;
  for (int i = 0; i < num_lists; ++i)
  {
    for (int j = 0; j < list_sizes[i]; ++j, ++k)
      for (int l = 0; l < width; ++l)
        merged_list[width*k+l] = lists_to_merge[i][width*j+l];
  }
  qsort(merged_list, k, width, compar);
}

static void parallel_qsort_with_regular_sampling(MPI_Comm comm, 
                                                 void* base, 
                                                 size_t nel, 
                                                 size_t width,
                                                 int (*compar)(const void* left, const void* right))
{
  // Size up our communicator.
  int nprocs, rank;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);

  // Now we begin by sorting local data.
  qsort(base, nel, width, compar);

  // On single process sorts, we're finished!
  if (nprocs == 1) return;

  // Here's our local array in bytes.
  uint8_t* base_bytes = base;

  // We will store sampled pivots in this buffer. NOTE: This is going to be 
  // the bottleneck for large numbers of processes!
  int pivot_buffer_size = nprocs * nprocs;
  uint8_t* pivot_buffer = polymec_malloc(width * pivot_buffer_size);

  // Now select pivots by sampling the sorted list.
  for (int p = 0; p < nprocs; ++p)
    for (int i = 0; i < width; ++i)
      pivot_buffer[p*width+i] = base_bytes[p*nel*width/nprocs+i];

  // The root process (rank 0) gathers all pivot candidates from the others.
  if (rank == 0)
    MPI_Gather(MPI_IN_PLACE, width*nprocs, MPI_BYTE, pivot_buffer, width*nprocs, MPI_BYTE, 0, comm);
  else
    MPI_Gather(pivot_buffer, width*nprocs, MPI_BYTE, pivot_buffer, width*nprocs, MPI_BYTE, 0, comm);

  // Now the root merges all of the various pivot buffers together and 
  // selects pivots to broadcast.
  if (rank == 0)
  {
    // Merge the pivot lists.
    uint8_t* pivot_lists[nprocs];
    int pivot_lengths[nprocs];
    for (int p = 0; p < nprocs; ++p)
    {
      pivot_lists[p] = &pivot_buffer[width*p];
      pivot_lengths[p] = nprocs;
    }
    uint8_t* merged_pivots = polymec_malloc(width * nprocs * nprocs);
    merge_sorted_lists(pivot_lists, pivot_lengths, nprocs, width, compar, merged_pivots);

    // Select pivots to broadcast as partition pivots for other processes.
    for (int p = 1; p < nprocs; ++p)
      for (int i = 0; i < width; ++i)
        pivot_buffer[(p-1)*width+i] = merged_pivots[p*nprocs*width+i];

    polymec_free(merged_pivots);
  }

  // The root process now broadcasts the pivots.
  MPI_Bcast(pivot_buffer, width*(nprocs-1), MPI_BYTE, 0, comm);

  // Now all processes divide their data into nprocs different classes. This 
  // will help us get the right data onto the right processes. 
  int class_offset[nprocs], class_size[nprocs], k = 0;
  for (int i = 0; i < nprocs-1; ++i)
  {
    class_offset[i] = k;
    class_size[i] = 0;

    // Populate this class with data.
    while ((k < nel) && (compar(&base_bytes[k*width], &pivot_buffer[i*width]) < 0))
    {
      ++(class_size[i]);
      ++k;
    }
  }
  class_offset[nprocs-1] = k;
  class_size[nprocs-1] = nel - k;

  // Now each process gathers the data belonging to its class.
  int data_sizes[nprocs], data_byte_sizes[nprocs];
  int data_byte_offsets[nprocs];
  uint8_t* class_data = polymec_malloc(width * nel);
  for (int p = 0; p < nprocs; ++p)
  {
    // Gather the various sizes of the data in this class from other 
    // processes to process p.
    MPI_Gather(&class_size[p], 1, MPI_INT, data_sizes, 1, MPI_INT, p, comm);
    for (int p = 0; p < nprocs; ++p)
      data_byte_sizes[p] = width * data_sizes[p];

    // Build a compressed representation of the offsets of the data
    // on process p.
    if (rank == p)
    {
      data_byte_offsets[0] = 0;
      for (int pp = 1; pp < nprocs; ++pp)
        data_byte_offsets[pp] = data_byte_offsets[pp-1] + data_byte_sizes[pp-1];
    }

    // Now process p gathers the data in this class.
    MPI_Gatherv(&base_bytes[width*class_offset[p]], class_size[p],
                MPI_BYTE, class_data, data_byte_sizes, data_byte_offsets, 
                MPI_BYTE, p, comm);
  }

  // Now each process merges the class data it has received.
  uint8_t* class_lists[nprocs];
  for (int p = 0; p < nprocs; ++p)
    class_lists[p] = &class_data[width * data_byte_offsets[p]];
  merge_sorted_lists(class_lists, data_sizes, nprocs, width, compar, base);

  // Clean up.
  polymec_free(pivot_buffer);
  polymec_free(class_data);
}

static void parallel_qsort_with_random_sampling(MPI_Comm comm, 
                                                void* base, 
                                                size_t nel, 
                                                size_t width,
                                                int (*compar)(const void* left, const void* right),
                                                rng_t* rng)
{
  POLYMEC_NOT_IMPLEMENTED;
}

void parallel_qsort_with_sampling(MPI_Comm comm, 
                                  void* base, 
                                  size_t nel, 
                                  size_t width,
                                  int (*compar)(const void* left, const void* right),
                                  rng_t* rng)
{
  if (rng == NULL)
    parallel_qsort_with_regular_sampling(comm, base, nel, width, compar);
  else
    parallel_qsort_with_random_sampling(comm, base, nel, width, compar, rng);
}

