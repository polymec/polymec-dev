// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmocka.h"
#include "core/exchanger.h"

static void test_migrator_new(void** state)
{
  migrator_t* m = migrator_new(MPI_COMM_WORLD);
  m = NULL;
}

static void test_migrator_from_global_partition(void** state)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  int nproc, rank;
  MPI_Comm_size(comm, &nproc);
  MPI_Comm_rank(comm, &rank);

  // Distribute 100 data to each process.
  int N = 100*nproc;
  int64_t P[N];
  for (int i = 0; i < N; ++i)
    P[i] = i / 100;
  migrator_t* m = migrator_from_global_partition(comm, P, N);

  // Migrate the data.
  real_t data[N];
  for (int i = 0; i < N; ++i)
    data[i] = 1.0 * i;
  int n = (rank == 0) ? N : 0;
  migrator_transfer(m, data, &n, 1, 0, MPI_REAL_T);
  assert_int_equal(n, 100);

  // Check it.
  for (int i = 0; i < n; ++i)
    assert_true(ABS(data[i] - (100.0 * rank + 1.0 * i)) < 1e-12);

  m = NULL;
}

static void test_migrator_from_local_partition_full_cycle(void** state)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  int nproc, rank;
  MPI_Comm_size(comm, &nproc);
  MPI_Comm_rank(comm, &rank);

  // Shift data from each process to the one on its "right."
  int N = 100;
  int64_t P[N];
  for (int i = 0; i < N; ++i)
    P[i] = (rank + 1) % nproc;
  migrator_t* m = migrator_from_local_partition(comm, P, N);

  // Migrate the data.
  real_t data[N];
  for (int i = 0; i < N; ++i)
    data[i] = 1.0 * (N * rank + i);
  int n = N;
  migrator_transfer(m, data, &n, 1, 0, MPI_REAL_T);
  assert_int_equal(n, N);

  // Check it.
  if (rank > 0)
  {
    for (int i = 0; i < n; ++i)
      assert_true(ABS(data[i] - 1.0 * (N * (rank-1) + i)) < 1e-12);
  }
  else
  {
    for (int i = 0; i < n; ++i)
      assert_true(ABS(data[i] - 1.0 * (N * (nproc-1) + i)) < 1e-12);
  }

  m = NULL;
}

static void test_migrator_from_local_partition_half_cycle(void** state)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  int nproc, rank;
  MPI_Comm_size(comm, &nproc);
  MPI_Comm_rank(comm, &rank);

  // Shift half of the data from each process to the one on its "right,"
  // and keep the remainder.
  int N = 100;
  int64_t P[N];
  for (int i = 0; i < N/2; ++i)
    P[i] = rank % nproc;
  for (int i = N/2; i < N; ++i)
    P[i] = (rank + 1) % nproc;
  migrator_t* m = migrator_from_local_partition(comm, P, N);

  // Migrate the data.
  real_t data[N];
  for (int i = 0; i < N; ++i)
    data[i] = 1.0 * rank;
  int n = N;
  migrator_transfer(m, data, &n, 1, 0, MPI_REAL_T);
  assert_int_equal(n, 100);

  m = NULL;
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_migrator_new),
    cmocka_unit_test(test_migrator_from_global_partition),
    cmocka_unit_test(test_migrator_from_local_partition_full_cycle),
    cmocka_unit_test(test_migrator_from_local_partition_half_cycle)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
