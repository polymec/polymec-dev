// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <stdlib.h>
#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmocka.h"
#include "core/rng.h"
#include "core/parallel_sort.h"

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

void test_sort_random_ints_with_fixed_size(void** state)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  rng_t* rng = host_rng_new();

  // Generate and sort a list of integers, 100 per process.
  size_t size = 100;
  int* ints = random_ints(rng, size);
  parallel_sort(comm, ints, size, sizeof(int), compare_ints);
  assert_true(list_is_ordered(comm, ints, size, sizeof(int), compare_ints));
  polymec_free(ints);
}

void test_sort_random_doubles_with_fixed_size(void** state)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  rng_t* rng = host_rng_new();

  // Generate and sort a list of doubles, 100 per process.
  size_t size = 100;
  double* doubles = random_doubles(rng, size);
  parallel_sort(comm, doubles, size, sizeof(double), compare_doubles);
  assert_true(list_is_ordered(comm, doubles, size, sizeof(double), compare_doubles));
  polymec_free(doubles);
}

void test_sort_random_ints_with_var_size(void** state)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  rng_t* rng = host_rng_new();

  // Generate and sort a list of integers, with variable sizes per process.
  // Even processes get 100, odd ones get 50.
  int rank;
  MPI_Comm_rank(comm, &rank);
  size_t size = (rank % 2 == 0) ? 100 : 50;
  int* ints = random_ints(rng, size);
  parallel_sort(comm, ints, size, sizeof(int), compare_ints);
  assert_true(list_is_ordered(comm, ints, size, sizeof(int), compare_ints));
  polymec_free(ints);
}

void test_sort_random_doubles_with_var_size(void** state)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  rng_t* rng = host_rng_new();

  // Generate and sort a list of doubles, with variable sizes per process.
  // Even processes get 100, odd ones get 50.
  int rank;
  MPI_Comm_rank(comm, &rank);
  size_t size = (rank % 2 == 0) ? 100 : 50;
  double* doubles = random_doubles(rng, size);
  parallel_sort(comm, doubles, size, sizeof(double), compare_doubles);
  assert_true(list_is_ordered(comm, doubles, size, sizeof(double), compare_doubles));
  polymec_free(doubles);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);

  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_sort_random_ints_with_fixed_size),
    cmocka_unit_test(test_sort_random_doubles_with_fixed_size),
    cmocka_unit_test(test_sort_random_ints_with_var_size),
    cmocka_unit_test(test_sort_random_doubles_with_var_size)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
