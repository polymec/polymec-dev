// Copyright (c) 2012-2019, Jeffrey N. Johnson
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

static void test_exchanger_new(void** state)
{
  exchanger_t* ex = exchanger_new(MPI_COMM_WORLD);
  release_ref(ex);
}

static void test_exchanger_construct(void** state)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  int nproc, rank;
  MPI_Comm_size(comm, &nproc);
  MPI_Comm_rank(comm, &rank);
  int N = 100*nproc;
  exchanger_t* ex = exchanger_new(comm);
  assert_true(exchanger_comm(ex) == comm);
  if (nproc > 1)
  {
    int send_indices[N/nproc];
    for (int i = 0; i < N/nproc; ++i)
      send_indices[i] = i;
    exchanger_set_send(ex, (rank+1) % nproc, send_indices, N/nproc, true);
    int receive_indices[N/nproc];
    for (int i = 0; i < N/nproc; ++i)
      receive_indices[i] = i;
    exchanger_set_receive(ex, (rank+nproc-1) % nproc, receive_indices, N/nproc, true);
  }
  release_ref(ex);
}

static void test_exchanger_construct_and_delete(void** state)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  int nproc, rank;
  MPI_Comm_size(comm, &nproc);
  MPI_Comm_rank(comm, &rank);
  int N = 100*nproc;
  exchanger_t* ex = exchanger_new(comm);
  exchanger_set_send_offset(ex, 0);
  exchanger_set_receive_offset(ex, 0);
  if (nproc > 1)
  {
    int send_indices[N/nproc];
    for (int i = 0; i < N/nproc; ++i)
      send_indices[i] = i;
    exchanger_set_send(ex, (rank+1) % nproc, send_indices, N/nproc, true);
    assert_true(exchanger_max_send(ex) != rank);
    int receive_indices[N/nproc];
    for (int i = 0; i < N/nproc; ++i)
      receive_indices[i] = i;
    exchanger_set_receive(ex, (rank+nproc-1) % nproc, receive_indices, N/nproc, true);
    assert_true(exchanger_max_receive(ex) != rank);

    exchanger_delete_send(ex, 1);
    exchanger_delete_receive(ex, 1);
  }
  release_ref(ex);
}

static void test_exchanger_is_valid_and_dl_detection(void** state)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  int nproc, rank;
  MPI_Comm_size(comm, &nproc);
  MPI_Comm_rank(comm, &rank);

  real_t grace_period = 0.2;

  // First we put together a bad exchanger and check that exchanger_is_valid
  // tells us it's bad, and that deadlock detection captures the deadlock.
  exchanger_t* ex = exchanger_new(comm);
  exchanger_enable_deadlock_detection(ex, grace_period, 0, stdout);
  assert_true(exchanger_deadlock_detection_enabled(ex));
  if (nproc > 1)
  {
    // We set up a receive from another rank, but we don't send anything.
    int receive_indices[100];
    for (int i = 0; i < 100; ++i)
      receive_indices[i] = i;
    exchanger_set_receive(ex, (rank+1) % nproc, receive_indices, 100, true);

    // Verify that it's a bad exchanger.
    char* reason;
    bool result = exchanger_is_valid(ex, &reason);
    assert_false(result);
    log_debug("Exchanger is invalid (expected): %s", reason);
    result = exchanger_is_valid(ex, NULL);
    assert_false(result);

    // Now try to exchange data.
    long data[100];
    for (int i = 0; i < 100; ++i)
      data[i] = (long)rank;
    exchanger_exchange(ex, data, 1, 0, MPI_LONG);
  }

  // Next, we put together a good exchanger and check that is_valid/deadlock
  // detection give it a pass.
  release_ref(ex);
  ex = exchanger_new(comm);
  exchanger_enable_deadlock_detection(ex, grace_period, 0, stdout);
  assert_true(exchanger_deadlock_detection_enabled(ex));
  if (nproc > 1)
  {
    // We set up a send/receive ring.
    int receive_indices[100];
    for (int i = 0; i < 100; ++i)
      receive_indices[i] = i;
    int receive_proc = (rank+1) % nproc;
    exchanger_set_receive(ex, receive_proc, receive_indices, 100, true);
    int send_indices[100];
    for (int i = 0; i < 100; ++i)
      send_indices[i] = i;
    int send_proc = ((rank-1) >= 0) ? rank-1 : rank+nproc-1;
    exchanger_set_send(ex, send_proc, send_indices, 100, true);

    // Verify that it's a good exchanger.
    char* reason;
    bool result = exchanger_is_valid(ex, &reason);
    assert_true(result);
    result = exchanger_is_valid(ex, NULL);
    assert_true(result);

    // Now try to exchange data.
#define EXCHANGE_DATA(c_type, mpi_type, expr) \
    { \
      c_type data[100]; \
      for (int i = 0; i < 100; ++i) \
        data[i] = (c_type)(expr); \
      exchanger_exchange(ex, data, 1, 0, mpi_type); \
      assert_true(data[0] != (c_type)rank); \
    }
    EXCHANGE_DATA(char, MPI_CHAR, rank)
    EXCHANGE_DATA(uint8_t, MPI_BYTE, rank)
    EXCHANGE_DATA(int, MPI_INT, rank)
    EXCHANGE_DATA(long, MPI_LONG, rank)
    EXCHANGE_DATA(long long, MPI_LONG_LONG, rank)
    EXCHANGE_DATA(uint64_t, MPI_UINT64_T, rank)
    EXCHANGE_DATA(int64_t, MPI_INT64_T, rank)
#undef EXCHANGE_DATA
  }
  release_ref(ex);
}

static void test_exchanger_local_copy(void** state)
{
  exchanger_t* ex = exchanger_new(MPI_COMM_SELF);

  // The exchanger copies values at even indices to odd ones, and values at odd
  // indices to even ones.
  int send_indices[100], receive_indices[100];
  for (int i = 0; i < 100; ++i)
  {
    send_indices[i] = i;
    if ((i % 2) != 0)
      receive_indices[i] = i-1;
    else
      receive_indices[i] = i+1;
  }
  exchanger_set_send(ex, 0, send_indices, 100, true);
  exchanger_set_receive(ex, 0, receive_indices, 100, true);

  // Verify that it's a good exchanger.
  bool result = exchanger_is_valid(ex, NULL);
  assert_true(result);

  // Now try to exchange data.
#define EXCHANGE_DATA(c_type, mpi_type, expr) \
  { \
    c_type data[100]; \
    for (int i = 0; i < 100; ++i) \
      data[i] = (c_type)i; \
    exchanger_exchange(ex, data, 1, 0, mpi_type); \
    for (int i = 0; i < 100; ++i) \
      assert_true(data[i] == ((i % 2) != 0) ? i-1 : i+1); \
  }
  EXCHANGE_DATA(char, MPI_CHAR, rank)
  EXCHANGE_DATA(uint8_t, MPI_BYTE, rank)
  EXCHANGE_DATA(int, MPI_INT, rank)
  EXCHANGE_DATA(long, MPI_LONG, rank)
  EXCHANGE_DATA(long long, MPI_LONG_LONG, rank)
  EXCHANGE_DATA(uint64_t, MPI_UINT64_T, rank)
  EXCHANGE_DATA(int64_t, MPI_INT64_T, rank)
#undef EXCHANGE_DATA
  release_ref(ex);
}

static void test_exchanger_reduce(void** state)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  int nproc, rank;
  MPI_Comm_size(comm, &nproc);
  MPI_Comm_rank(comm, &rank);

  exchanger_t* ex = exchanger_new(comm);

  // All processes have one index (0), and all send their values to each
  // other for various reductions.
  int index = 0;
  for (int proc = 0; proc < nproc; ++proc)
    exchanger_set_receive(ex, proc, &index, 1, true);
  for (int proc = 0; proc < nproc; ++proc)
    exchanger_set_send(ex, proc, &index, 1, true);
  assert_true(exchanger_aggregates_data(ex) || (nproc == 1));

  // Verify that it's a good exchanger.
  bool result = exchanger_is_valid(ex, NULL);
  assert_true(result);

  // Now try to exchange data.
#define EXCHANGE_AND_REDUCE_DATA(c_type, mpi_type, reduce_op, start, answer) \
  { \
    c_type data[1]; \
    data[0] = (c_type)(start); \
    exchanger_set_reducer(ex, reduce_op); \
    exchanger_exchange(ex, data, 1, 0, mpi_type); \
    assert_true(data[0] == (c_type)answer); \
  }

#define TEST_ALL_REDUCERS(c_type, mpi_type) \
  { \
    c_type start = (c_type)(rank+1); \
    c_type sum = 0; \
    c_type product = (c_type)1; \
    for (int p = 0; p < nproc; ++p) \
    { \
      sum += (p+1); \
      product *= (p+1); \
    } \
    c_type min = (c_type)1; \
    c_type min_rank = (c_type)1; \
    c_type max = (c_type)nproc; \
    c_type max_rank = (c_type)nproc; \
    EXCHANGE_AND_REDUCE_DATA(c_type, mpi_type, EXCHANGER_SUM, start, sum) \
    EXCHANGE_AND_REDUCE_DATA(c_type, mpi_type, EXCHANGER_PRODUCT, start, product) \
    EXCHANGE_AND_REDUCE_DATA(c_type, mpi_type, EXCHANGER_MIN, start, min) \
    EXCHANGE_AND_REDUCE_DATA(c_type, mpi_type, EXCHANGER_MAX, start, max) \
    EXCHANGE_AND_REDUCE_DATA(c_type, mpi_type, EXCHANGER_MIN_RANK, start, min_rank) \
    EXCHANGE_AND_REDUCE_DATA(c_type, mpi_type, EXCHANGER_MAX_RANK, start, max_rank) \
  }

  TEST_ALL_REDUCERS(char, MPI_CHAR)
  TEST_ALL_REDUCERS(uint8_t, MPI_BYTE)
  TEST_ALL_REDUCERS(int, MPI_INT)
  TEST_ALL_REDUCERS(long, MPI_LONG)
  TEST_ALL_REDUCERS(long long, MPI_LONG_LONG)
  TEST_ALL_REDUCERS(uint64_t, MPI_UINT64_T)
  TEST_ALL_REDUCERS(int64_t, MPI_INT64_T)
#undef EXCHANGE_AND_REDUCE_DATA
#undef TEST_ALL_REDUCERS
  release_ref(ex);
}

int main(int argc, char* argv[])
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] =
  {
    cmocka_unit_test(test_exchanger_new),
    cmocka_unit_test(test_exchanger_construct),
    cmocka_unit_test(test_exchanger_construct_and_delete),
    cmocka_unit_test(test_exchanger_is_valid_and_dl_detection),
    cmocka_unit_test(test_exchanger_local_copy),
    cmocka_unit_test(test_exchanger_reduce)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
