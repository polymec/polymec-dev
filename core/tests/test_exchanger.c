// Copyright (c) 2012-2017, Jeffrey N. Johnson
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
  exchanger_t* exchanger = exchanger_new(MPI_COMM_WORLD);
  exchanger = NULL;
}

static void test_exchanger_construct(void** state)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  int nproc, rank;
  MPI_Comm_size(comm, &nproc);
  MPI_Comm_rank(comm, &rank);
  int N = 100*nproc;
  exchanger_t* exchanger = exchanger_new(comm);
  assert_true(exchanger_comm(exchanger) == comm);
  if (nproc > 1)
  {
    int send_indices[N/nproc];
    for (int i = 0; i < N/nproc; ++i)
      send_indices[i] = i;
    exchanger_set_send(exchanger, (rank+1) % nproc, send_indices, N/nproc, true);
    int receive_indices[N/nproc];
    for (int i = 0; i < N/nproc; ++i)
      receive_indices[i] = i;
    exchanger_set_receive(exchanger, (rank+nproc-1) % nproc, receive_indices, N/nproc, true);
  }
  exchanger = NULL;
}

static void test_exchanger_construct_and_delete(void** state)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  int nproc, rank;
  MPI_Comm_size(comm, &nproc);
  MPI_Comm_rank(comm, &rank);
  int N = 100*nproc;
  exchanger_t* exchanger = exchanger_new(comm);
  exchanger_set_send_offset(exchanger, 0);
  exchanger_set_receive_offset(exchanger, 0);
  if (nproc > 1)
  {
    int send_indices[N/nproc];
    for (int i = 0; i < N/nproc; ++i)
      send_indices[i] = i;
    exchanger_set_send(exchanger, (rank+1) % nproc, send_indices, N/nproc, true);
    assert_true(exchanger_max_send(exchanger) != rank);
    int receive_indices[N/nproc];
    for (int i = 0; i < N/nproc; ++i)
      receive_indices[i] = i;
    exchanger_set_receive(exchanger, (rank+nproc-1) % nproc, receive_indices, N/nproc, true);
    assert_true(exchanger_max_receive(exchanger) != rank);

    exchanger_delete_send(exchanger, 1);
    exchanger_delete_receive(exchanger, 1);
  }
  exchanger = NULL;
}

static bool _bad_exchanger = false;

static void handle_bad_exchanger(const char* format, ...)
{
  printf("That there exchanger you got is broked!");
  _bad_exchanger = true;
}

static void test_exchanger_verify_and_dl_detection(void** state)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  int nproc, rank;
  MPI_Comm_size(comm, &nproc);
  MPI_Comm_rank(comm, &rank);

  real_t grace_period = 0.2;

  // First we put together a bad exchanger and check that verify tells us 
  // it's bad, and that deadlock detection captures the deadlock.
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
    exchanger_verify(ex, handle_bad_exchanger);
//    assert_true(_bad_exchanger); // FIXME
    _bad_exchanger = false;

    // Now try to exchange data.
    long data[100];
    for (int i = 0; i < 100; ++i)
      data[i] = (long)rank;
    exchanger_exchange(ex, data, 1, 0, MPI_LONG);
//printf("%d %d\n", (int)data[0], rank);
//    assert_true(data[0] == (long)rank); // FIXME
  }

  // Next, we put together a good exchanger and check that verify/deadlock 
  // detection give it a pass.
  ex = exchanger_new(comm);
  exchanger_enable_deadlock_detection(ex, grace_period, 0, stdout);
  assert_true(exchanger_deadlock_detection_enabled(ex));
  if (nproc > 1)
  {
    // We set up a send/receive ring.
    int receive_indices[100];
    for (int i = 0; i < 100; ++i)
      receive_indices[i] = i;
    exchanger_set_receive(ex, (rank+1) % nproc, receive_indices, 100, true);
    int send_indices[100];
    for (int i = 0; i < 100; ++i)
      send_indices[i] = i;
    int send_proc = ((rank-1) >= 0) ? rank-1 : rank+nproc-1;
    exchanger_set_send(ex, send_proc, send_indices, 100, true);

    // Verify that it's a good exchanger.
    exchanger_verify(ex, handle_bad_exchanger);
    assert_false(_bad_exchanger);

    // Now try to exchange data.
    long long data[100];
    for (int i = 0; i < 100; ++i)
      data[i] = (long long)rank;
    exchanger_exchange(ex, data, 1, 0, MPI_LONG_LONG);
    assert_true(data[0] != (long long)rank);
  }
  ex = NULL;

}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_exchanger_new),
    cmocka_unit_test(test_exchanger_construct),
    cmocka_unit_test(test_exchanger_construct_and_delete),
    cmocka_unit_test(test_exchanger_verify_and_dl_detection)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
