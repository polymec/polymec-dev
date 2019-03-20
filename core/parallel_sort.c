// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/parallel_sort.h"

#if POLYMEC_HAVE_MPI

// Returns the index of the maximum element in the given array.
static int max_index(void* data,
                     index_t nel,
                     size_t width,
                     int (*comp)(const void* left, const void* right))
{
  void* max = data;
  int i_max = 0;
  int8_t* bytes = data;

  for (int i = 1; i < nel; ++i)
  {
    if (comp(&bytes[i*width], max) > 0)
    {
      max = &bytes[i*width];
      i_max = i;
    }
  }
  return i_max;
}

// Returns the index of the minimum element in the given array.
static int min_index(void* data,
                     index_t nel,
                     size_t width,
                     int (*comp)(const void* left, const void* right))
{
  void* min = data;
  int i_min = 0;
  int8_t* bytes = data;

  for (int i = 1; i < nel; ++i)
  {
    if (comp(&bytes[i*width], min) < 0)
    {
      min = &bytes[i*width];
      i_min = i;
    }
  }
  return i_min;
}

// Swaps two sets of bytes of the given size.
static void swap_bytes(void* x, void* y, size_t size)
{
  int8_t temp[size];
  memcpy(temp, x, size);
  memcpy(x, y, size);
  memcpy(y, temp, size);
}

#endif

// This parallel sort is adapted from
// http://cs.umw.edu/~finlayson/class/fall14/cpsc425/notes/18-sorting.html.
void parallel_sort(MPI_Comm comm,
                   void* base,
                   size_t nel,
                   size_t width,
                   int (*comp)(const void* left, const void* right))
{
#if POLYMEC_HAVE_MPI
  int rank, nprocs;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &nprocs);

  // If we're only on 1 proc, we can just do a local sort.
  if (nprocs == 1)
  {
    qsort(base, nel, width, comp);
    return;
  }

  // Find our partners for even and odd phases.
  int even_partner = ((rank % 2) == 0) ? rank + 1 : rank - 1;
  int odd_partner = ((rank % 2) == 0) ? rank - 1 : rank + 1;

  // Get the sizes for our even- and odd-phase partner processes.
  index_t N = (index_t)nel;
  index_t N_even = 0, N_odd = 0;
  if ((rank % 2) == 0)
  {
    if ((even_partner >= 0) && (even_partner < nprocs))
    {
      MPI_Send(&N, 1, MPI_INDEX_T, even_partner, 0, comm);
      MPI_Recv(&N_even, 1, MPI_INDEX_T, even_partner, 0, comm, MPI_STATUS_IGNORE);
    }
    if ((odd_partner >= 0) && (odd_partner < nprocs))
    {
      MPI_Recv(&N_odd, 1, MPI_INDEX_T, odd_partner, 0, comm, MPI_STATUS_IGNORE);
      MPI_Send(&N, 1, MPI_INDEX_T, odd_partner, 0, comm);
    }
  }
  else
  {
    if ((even_partner >= 0) && (even_partner < nprocs))
    {
      MPI_Recv(&N_even, 1, MPI_INDEX_T, even_partner, 0, comm, MPI_STATUS_IGNORE);
      MPI_Send(&N, 1, MPI_INDEX_T, even_partner, 0, comm);
    }
    if ((odd_partner >= 0) && (odd_partner < nprocs))
    {
      MPI_Send(&N, 1, MPI_INDEX_T, odd_partner, 0, comm);
      MPI_Recv(&N_odd, 1, MPI_INDEX_T, odd_partner, 0, comm, MPI_STATUS_IGNORE);
    }
  }

  // Compute the ratio of the max to min numbers of elements for all process.
  // The number of phases in the sort is this ratio times nprocs.
  int N_max, N_min;
  MPI_Allreduce(&N, &N_max, 1, MPI_INT, MPI_MAX, comm);
  MPI_Allreduce(&N, &N_min, 1, MPI_INT, MPI_MIN, comm);
  real_t max_min_ratio = 1.0 * N_max / N_min;
  int nphases = (int)(ceil(max_min_ratio * nprocs));

  // Allocate storage for the exchange.
  void* other = polymec_malloc(width * MAX(N_even, N_odd));

  // Sort in nprocs phases (even and odd, depending upon i's parity).
  for (int i = 0; i < nphases; ++i)
  {
    // Sort our local data.
    qsort(base, N, width, comp);

    // Find our sorting "partner" for this phase and its size (N_partner).
    int N_partner, partner;
    if (i % 2 == 0) // even phase
    {
      N_partner = (int)N_even;
      partner = even_partner;
    }
    else // odd phase
    {
      N_partner = (int)N_odd;
      partner = odd_partner;
    }

    // If our partner process falls off one of the ends, skip this phase.
    if ((partner < 0) || (partner >= nprocs))
      continue;

    // Exchange data - even processes send first and odd processes receive
    // first, avoiding deadlock conditions.
    if (rank % 2 == 0)
    {
      MPI_Send(base, (int)(N*width), MPI_BYTE, partner, 0, MPI_COMM_WORLD);
      MPI_Recv(other, (int)(N_partner*width), MPI_BYTE, partner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    else
    {
      MPI_Recv(other, (int)(N_partner*width), MPI_BYTE, partner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Send(base, (int)(N*width), MPI_BYTE, partner, 0, MPI_COMM_WORLD);
    }

    // Merge base and other based on our preferred range of values.
    if (rank < partner)
    {
      // Seek smaller values.
      int num_swaps = 0;
      while (num_swaps < MIN(N, N_partner))
      {
        // Find the index of the smallest element in the other array.
        int i_min = min_index(other, N_partner, width, comp);

        // Find the index of the largest element in our array.
        int i_max = max_index(base, N, width, comp);

        // If the smallest one in the other array is less than the largest
        // in ours, swap them.
        int8_t* base_bytes = base;
        int8_t* other_bytes = other;
        if (comp(&other_bytes[width*i_min], &base_bytes[width*i_max]) < 0)
        {
          swap_bytes(&other_bytes[width*i_min], &base_bytes[width*i_max], width);
          ++num_swaps;
        }
        else
        {
          // The smallest are now in base.
          break;
        }
      }
    }
    else
    {
      // Seek larger values.
      int num_swaps = 0;
      while (num_swaps < MIN(N, N_partner))
      {
        // Find the index of the largest element in the other array.
        int i_max = max_index(other, N_partner, width, comp);

        // Find the index of the smallest element in our array.
        int i_min = min_index(base, N, width, comp);

        // If the largest one in the other array is bigger than the smallest
        // in ours, swap them.
        int8_t* base_bytes = base;
        int8_t* other_bytes = other;
        if (comp(&other_bytes[width*i_max], &base_bytes[width*i_min]) > 0)
        {
          swap_bytes(&other_bytes[width*i_max], &base_bytes[width*i_min], width);
          ++num_swaps;
        }
        else
        {
          // The largest are now in base.
          break;
        }
      }
    }
  }

  // Clean up.
  polymec_free(other);

  // Sort one last time.
  qsort(base, nel, width, comp);

#else
  // Just sort our local data.
  qsort(base, nel, width, comp);
#endif
}

