// Copyright (c) 2012-2015, Jeffrey N. Johnson
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

#include <stdlib.h>
#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmockery.h"
#include "core/parallel_qsort.h"

static int compare_ints(const void* l, const void* r)
{
  const int* left = l;
  const int* right = r;
  return (left[0] < right[0]) ? -1 : (left[0] == right[0]) ? 0 : 1;
}

static int compare_doubles(const void* l, const void* r)
{
  const double* left = l;
  const double* right = r;
  return (left[0] < right[0]) ? -1 : (left[0] == right[0]) ? 0 : 1;
}

static int* random_ints(rng_t* rng, size_t size)
{
  int* list = polymec_malloc(sizeof(int) * size);
  for (int i = 0; i < size; ++i)
    list[i] = rng_uniform_int(rng, 10000);
  return list;
}

static double* random_doubles(rng_t* rng, size_t size)
{
  double* list = polymec_malloc(sizeof(double) * size);
  for (int i = 0; i < size; ++i)
    list[i] = rng_uniform(rng);
  return list;
}

static bool list_is_ordered(MPI_Comm comm, 
                            void* base, 
                            size_t nel, 
                            size_t width, 
                            int (*compar)(const void* left, const void* right))
{
  int nprocs, rank;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);

  // Exchange endpoints with our neighbors.
  char* list = base;
  char left[width], right[width];
  MPI_Status statuses[2];
  MPI_Request requests[2];
  int num_requests = 0;
  if (rank > 0)
  {
    MPI_Irecv(left, width, MPI_BYTE, rank-1, 0, comm, &requests[num_requests]);
    ++num_requests;
  }
  if (rank < (nprocs - 1))
  {
    MPI_Irecv(right, width, MPI_BYTE, rank+1, 0, comm, &requests[num_requests]);
    ++num_requests;
  }
  if (rank > 0)
    MPI_Send(list, width, MPI_BYTE, rank-1, 0, comm);
  if (rank < (nprocs - 1))
    MPI_Send(&list[width*(nel-1)], width, MPI_BYTE, rank+1, 0, comm);

  MPI_Waitall(num_requests, requests, statuses);

  // Now make sure everything is properly sorted.
  if ((rank > 0) && (compar(left, &list[0]) > 0))
  {
    printf("list_is_ordered: last element on rank %d exceeds first on rank %d\n", rank-1, rank);
    return false;
  }

  for (int i = 1; i < nel; ++i)
  {
    if (compar(&list[width*(i-1)], &list[width*i]) > 0)
    {
      printf("list_is_ordered: element %d exceeds element %d on rank %d\n", i-1, i, rank);
      return false;
    }
  }

  if ((rank < (nprocs - 1)) && (compar(&list[width*(nel-1)], right) > 0))
  {
    printf("list_is_ordered: last element on rank %d exceeds first on rank %d\n", rank, rank+1);
    return false;
  }

  return true;
}

void test_regular_sampling(void** state)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  rng_t* rng = host_rng_new();

  // Generate and sort a list of integers, 100 per process.
  {
    size_t size = 100;
    int* ints = random_ints(rng, size);
int rank;
MPI_Comm_rank(comm, &rank);
//printf("%d: ", rank);
//for (int i = 0; i < size; ++i)
//  printf("%d ", ints[i]);
//printf("\n");
    parallel_qsort(comm, ints, size, sizeof(int), compare_ints, NULL);
printf("%d (%d): ", rank, size);
for (int i = 0; i < size; ++i)
  printf("%d ", ints[i]);
printf("\n");
    assert_true(list_is_ordered(comm, ints, size, sizeof(int), compare_ints));
    polymec_free(ints);
  }

  // Generate and sort a list of doubles, 100 per process.
return;
  {
    size_t size = 100;
    double* doubles = random_doubles(rng, size);
    parallel_qsort(comm, doubles, size, sizeof(double), compare_doubles, NULL);
    assert_true(list_is_ordered(comm, doubles, size, sizeof(double), compare_doubles));
    polymec_free(doubles);
  }
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);

  const UnitTest tests[] = 
  {
    unit_test(test_regular_sampling)
  };
  return run_tests(tests);
}
