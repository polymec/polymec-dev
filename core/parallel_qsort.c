// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/parallel_qsort.h"

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

static void reshuffle_if_necessary(MPI_Comm comm, 
                                   size_t width, 
                                   size_t new_size, 
                                   void* temp,
                                   size_t nel,
                                   void* base)
{
  int nprocs, rank;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);

  // If we have disturbed the number of elements per process, we reshuffle 
  // them now. This is potentially expensive, but much less confusing than 
  // imbalancing the load on each process for our purposes. Everyone needs 
  // to participate at this stage of the game, in any case.
  int reshuffle = (new_size == nel) ? 0 : 1, global_reshuffle;
  MPI_Allreduce(&reshuffle, &global_reshuffle, 1, MPI_INT, MPI_MAX, comm);
  if (global_reshuffle)
  {
    uint8_t* temp_bytes = temp;
    int num_excess_elems = new_size - nel;
    int num_elems_from_left = 0, num_elems_from_right = 0, num_elems_to_left = 0;
    uint8_t* elems_from_left = NULL;
    uint8_t* elems_from_right = NULL;

    // Get any elements we need from our neighbor.
    MPI_Status status;

    // Interact with the left neighbor.
printf("%d has %d excess elements.\n", rank, num_excess_elems);
    if (rank > 0)
    {
      MPI_Recv(&num_elems_from_left, 1, MPI_INT, rank-1, 0, comm, &status);
printf("%d is getting %d elems from the left.\n", rank, num_elems_from_left);
      num_excess_elems += num_elems_from_left;
      if (num_elems_from_left > 0)
      {
        elems_from_left = polymec_malloc(width * num_elems_from_left);
        MPI_Recv(elems_from_left, width * num_elems_from_left, MPI_BYTE, rank-1, 0, comm, &status);
      }
      else if (num_elems_from_left < 0)
      {
        num_elems_to_left = -num_elems_from_left;
        uint8_t* elems_to_left = polymec_malloc(width * num_elems_to_left);
        memcpy(elems_to_left, temp_bytes, width * num_elems_to_left);
        MPI_Send(elems_to_left, width * num_elems_to_left, MPI_BYTE, rank-1, 0, comm);
        polymec_free(elems_to_left);
      }
    }
printf("%d now has %d excess elements.\n", rank, num_excess_elems);

    // Interact with the right neighbor.
    if (rank < (nprocs - 1))
    {
      MPI_Send(&num_excess_elems, 1, MPI_INT, rank+1, 0, comm);
      if (num_excess_elems < 0)
      {
        // We are getting elements from our right neighbor.
        num_elems_from_right = -num_excess_elems;
printf("%d is getting %d elems from the right.\n", rank, num_elems_from_right);
        elems_from_right = polymec_malloc(width * num_elems_from_right);
        MPI_Recv(elems_from_right, width * num_elems_from_right, MPI_BYTE, rank+1, 0, comm, &status);
      }
      else if (num_excess_elems > 0)
      {
        // We are giving elements to our right neighbor.
        uint8_t* elems_to_right = polymec_malloc(width * num_excess_elems);
        if ((num_elems_from_left > 0) && (num_elems_from_left > nel))
        {
          // We are passing along excess elements we got from the left, since we got so many.
          int num_passed_on_from_left = num_elems_from_left - nel;
          memcpy(elems_to_right, elems_from_left[width*nel], width * num_passed_on_from_left);
          if (num_excess_elems > num_passed_on_from_left)
            memcpy(&elems_to_right[width*num_passed_on_from_left], &temp_bytes[num_elems_to_left], width * (num_excess_elems - num_passed_on_from_left));
          num_elems_from_left = nel;
        }
        else
          memcpy(elems_to_right, &temp_bytes[width*(nel+num_elems_to_left)], width * num_excess_elems);
        MPI_Send(elems_to_right, width * num_excess_elems, MPI_BYTE, rank+1, 0, comm);
        polymec_free(elems_to_right);
      }
    }

    // Copy the elements into place.
    uint8_t* base_bytes = base;
    memcpy(base_bytes, temp_bytes, width * nel);
    if (num_elems_from_left > 0)
    {
printf("%d: Getting stuff from left!\n", rank);
      if (num_elems_from_left > 0)
      {
        if (num_elems_from_left < nel)
          memcpy(&base_bytes[width * num_elems_from_left], temp_bytes, width * (nel - num_elems_from_left));
        memcpy(base_bytes, elems_from_left, width * MIN(nel, num_elems_from_left));
      }
    }
    else if (num_elems_to_left > 0)
    {
      if (num_elems_to_left < nel)
        memcpy(base_bytes, &temp_bytes[width * num_elems_to_left], width * nel);
    }
    if (num_elems_from_right > 0)
    {
printf("%d: Getting stuff from right!\n", rank);
      ASSERT(num_elems_from_right <= nel);
      memcpy(&base_bytes[width * (nel-num_elems_from_right)], elems_from_right, width * num_elems_from_right);
    }
  }
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
//{printf("%d: Selected pivots: [", rank);
//int* a = (int*)pivot_buffer;
//for (int i = 0; i < nprocs; ++i)
//printf("%d ", a[i]);
//printf("]\n");}

  // The root process (rank 0) gathers all pivot candidates from the others.
  if (rank == 0)
    MPI_Gather(MPI_IN_PLACE, width*nprocs, MPI_BYTE, pivot_buffer, width*nprocs, MPI_BYTE, 0, comm);
  else
    MPI_Gather(pivot_buffer, width*nprocs, MPI_BYTE, pivot_buffer, width*nprocs, MPI_BYTE, 0, comm);

  // Now the root merges all of the various pivot buffers together and 
  // selects pivots to broadcast.
  if (rank == 0)
  {
//{printf("Gathered pivots: [");
//int* a = (int*)pivot_buffer;
//for (int i = 0; i < nprocs*nprocs; ++i)
//printf("%d ", a[i]);
//printf("]\n");}
    // Merge the pivot lists.
    uint8_t* pivot_lists[nprocs];
    int pivot_lengths[nprocs];
    for (int p = 0; p < nprocs; ++p)
    {
      pivot_lists[p] = &pivot_buffer[width*p*nprocs];
      pivot_lengths[p] = nprocs;
    }
    uint8_t* merged_pivots = polymec_malloc(width * nprocs * nprocs);
    merge_sorted_lists(pivot_lists, pivot_lengths, nprocs, width, compar, merged_pivots);

    // Select pivots to broadcast as partition pivots for other processes.
    for (int p = 1; p < nprocs; ++p)
      for (int i = 0; i < width; ++i)
        pivot_buffer[(p-1)*width+i] = merged_pivots[p*nprocs*width+i];
//{printf("Merged pivots: [");
//int* a = (int*)pivot_buffer;
//for (int i = 0; i < nprocs; ++i)
//printf("%d ", a[i]);
//printf("]\n");}

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
    while ((k < nel) && (compar(&base_bytes[k*width], &pivot_buffer[i*width]) <= 0))
    {
      ++(class_size[i]);
      ++k;
    }
//printf("%d: %d members of class %d (starting at %d)\n", rank, class_size[i], i, class_offset[i]);
//{
//int* a = (int*)base_bytes;
//printf("%d: [", rank);
//for (int j = 0; j < class_size[i]; ++j)
//printf("%d ", a[class_offset[i]+j]);
//printf("]\n");
//}
  }
  class_offset[nprocs-1] = k;
  class_size[nprocs-1] = nel - k;
//printf("%d: %d members of class %d (starting at %d)\n", rank, class_size[nprocs-1], nprocs-1, class_offset[nprocs-1]);
//{
//int* a = (int*)base_bytes;
//printf("%d: [", rank);
//for (int j = 0; j < class_size[nprocs-1]; ++j)
//printf("%d ", a[class_offset[nprocs-1]+j]);
//printf("]\n");
//}

  // Gather the various sizes of the data in this class from other 
  // processes to process p.
  int data_sizes[nprocs], data_byte_sizes[nprocs];
  int data_byte_offsets[nprocs];
  for (int p = 0; p < nprocs; ++p)
  {
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
  }
  size_t new_size = (data_byte_offsets[nprocs-1] + data_byte_sizes[nprocs-1])/width;

  // Now each process gathers the data belonging to its class.
  uint8_t* class_data = polymec_malloc(width * new_size);
  for (int p = 0; p < nprocs; ++p)
  {
    // Now process p gathers the data in this class.
    MPI_Gatherv(&base_bytes[width*class_offset[p]], width * class_size[p],
                MPI_BYTE, class_data, data_byte_sizes, data_byte_offsets, 
                MPI_BYTE, p, comm);
  }
//{
//int* a = (int*)class_data;
//printf("%d: Unmerged class data: [", rank);
//for (int i = 0; i < nel; ++i)
//printf("%d ", a[i]);
//printf("]\n");
//}

  // If the new local size is not the old size, create a temporary local 
  // array.
  void* temp = base;
  if (new_size != nel)
    temp = polymec_malloc(width * new_size);

  // Now each process merges the class data it has received.
  uint8_t* class_lists[nprocs];
  for (int p = 0; p < nprocs; ++p)
    class_lists[p] = &class_data[data_byte_offsets[p]];
  merge_sorted_lists(class_lists, data_sizes, nprocs, width, compar, temp);
{
int* a = (int*)temp;
printf("%d: Unshuffled list: [", rank);
for (int i = 0; i < new_size; ++i)
printf("%d ", a[i]);
printf("]\n");
}

  reshuffle_if_necessary(comm, width, new_size, temp, nel, base);

  // Clean up.
  if (temp != base)
    polymec_free(temp);
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

void parallel_qsort(MPI_Comm comm, 
                    void* base, 
                    size_t nel, 
                    size_t width,
                    int (*compar)(const void* left, const void* right),
                    rng_t* rng)
{
  ASSERT(base != NULL);

  if (rng == NULL)
    parallel_qsort_with_regular_sampling(comm, base, nel, width, compar);
  else
    parallel_qsort_with_random_sampling(comm, base, nel, width, compar, rng);
}

