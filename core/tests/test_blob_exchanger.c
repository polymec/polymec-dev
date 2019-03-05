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
#include "core/blob_exchanger.h"

// Blobs of various types.
typedef struct
{
  char name[16];
  int age;
} small_widget_t;

typedef struct
{
  char name[16];
  int age;
  real_t spectrum[100];
} big_widget_t;

// This sets up processes in a ring for running tests.
static blob_exchanger_t* ring_exchanger(void** state)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  int nproc, rank;
  MPI_Comm_size(comm, &nproc);
  MPI_Comm_rank(comm, &rank);

  // We connect our processes in a ring. Everyone passes one small widget to
  // its left and one big widget to its right. The small widget has index 0
  // and the big one has index 1.
  blob_exchanger_proc_map_t* send_map = blob_exchanger_proc_map_new();
  int left = (rank > 0) ? rank - 1 : nproc - 1;
  int right = (rank < (nproc-1)) ? rank + 1 : 0;
  int small_blob = 0, big_blob = 1;

  blob_exchanger_proc_map_add_index(send_map, left, small_blob);
  blob_exchanger_proc_map_add_index(send_map, right, big_blob);

  // And of course we expect to receive blobs from our neighbors.
  blob_exchanger_proc_map_t* recv_map = blob_exchanger_proc_map_new();
  blob_exchanger_proc_map_add_index(recv_map, right, small_blob);
  blob_exchanger_proc_map_add_index(recv_map, left, big_blob);

  // Now we assign sizes to the blobs.
  blob_exchanger_size_map_t* size_map = blob_exchanger_size_map_new();
  blob_exchanger_size_map_insert(size_map, small_blob, sizeof(small_widget_t));
  blob_exchanger_size_map_insert(size_map, big_blob, sizeof(big_widget_t));

  // Create the exchanger.
  return blob_exchanger_new(comm, send_map, recv_map, size_map);
}

static void test_blob_exchanger_construct(void** state)
{
  blob_exchanger_t* ex = ring_exchanger(state);
  assert_true(blob_exchanger_comm(ex) == MPI_COMM_WORLD);

  // Create a blob buffer that can be used with this exchanger.
  blob_buffer_t* buffer = blob_exchanger_create_buffer(ex, 1);

  // Clean up.
  blob_buffer_free(buffer);
  release_ref(ex);
}

static void test_blob_exchanger_exchange(void** state)
{
  // Make our ring exchanger.
  blob_exchanger_t* ex = ring_exchanger(state);
  blob_exchanger_fprintf(ex, stdout);

  // Create a blob buffer that can be used with this exchanger.
  blob_buffer_t* buffer = blob_exchanger_create_buffer(ex, 1);

  // Initialize small and big widgets to exchange and copy them into
  // our buffer.
  int small_blob = 0, big_blob = 1;
  small_widget_t smallw;
  strcpy(smallw.name, "lil widget!");
  smallw.age = 10;
  blob_exchanger_copy_in(ex, small_blob, 1, &smallw, buffer);
  big_widget_t bigw;
  strcpy(bigw.name, "BIG widget!");
  bigw.age = 100;
  for (int i = 0; i < 100; ++i)
    bigw.spectrum[i] = 1.0 * i;
  blob_exchanger_copy_in(ex, big_blob, 1, &bigw, buffer);

  // Do the exchange.
  blob_exchanger_exchange(ex, 0, buffer);

  // Now get our received data.
  small_widget_t smallw1;
  blob_exchanger_copy_out(ex, buffer, small_blob, 1, &smallw1);
  big_widget_t bigw1;
  blob_exchanger_copy_out(ex, buffer, big_blob, 1, &bigw1);

  // Did we get what we expected?
  assert_int_equal(0, strcmp(smallw.name, "lil widget!"));
  assert_int_equal(10, smallw.age);
  assert_int_equal(0, strcmp(bigw.name, "BIG widget!"));
  assert_int_equal(100, bigw.age);
  for (int i = 0; i < 100; ++i)
    assert_true(reals_equal(1.0*i, bigw.spectrum[i]));

  // Clean up.
  blob_buffer_free(buffer);
  release_ref(ex);
}

/*
// This sets up a bad exchanger to use for deadlock detection.
static blob_exchanger_t* bad_exchanger(void** state)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  int nproc, rank;
  MPI_Comm_size(comm, &nproc);
  MPI_Comm_rank(comm, &rank);

  // As in our ring exchanger, everyone passes one small widget to its left and
  // one big widget to its right. The small widget has index 0 and the big one
  // has index 1.
  blob_exchanger_proc_map_t* send_map = blob_exchanger_proc_map_new();
  int left = (rank > 0) ? rank - 1 : nproc - 1;
  int right = (rank < (nproc-1)) ? rank + 1 : 0;
  int small_blob = 0, big_blob = 1;
  blob_exchanger_proc_map_add_index(send_map, left, small_blob);
  blob_exchanger_proc_map_add_index(send_map, right, big_blob);

  // Here's a flawed receive map.
  blob_exchanger_proc_map_t* recv_map = blob_exchanger_proc_map_new();
  blob_exchanger_proc_map_add_index(recv_map, right, big_blob);
  blob_exchanger_proc_map_add_index(recv_map, right, big_blob+1);
  blob_exchanger_proc_map_add_index(recv_map, right, big_blob+2);

  // Assign blob sizes. La dee da... everything's *fine*!
  blob_exchanger_size_map_t* size_map = blob_exchanger_size_map_new();
  blob_exchanger_size_map_insert(size_map, small_blob, sizeof(small_widget_t));
  blob_exchanger_size_map_insert(size_map, big_blob, sizeof(big_widget_t));
  blob_exchanger_size_map_insert(size_map, big_blob+1, sizeof(big_widget_t));
  blob_exchanger_size_map_insert(size_map, big_blob+2, sizeof(big_widget_t));

  // Create the exchanger.
  return blob_exchanger_new(comm, send_map, recv_map, size_map);
}

static void test_blob_exchanger_verify_and_dl_detection(void** state)
{
  // We skip this for single process runs.
  int nproc;
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  if (nproc == 1)
    return;

  // Create a busted exchanger.
  blob_exchanger_t* ex = bad_exchanger(state);

  // Enable deadlock detection.
  real_t grace_period = 0.2;
  blob_exchanger_enable_deadlock_detection(ex, grace_period, 0, stdout);
  assert_true(blob_exchanger_deadlock_detection_enabled(ex));

  // Verify that it's a bad exchanger.
  bool result = blob_exchanger_verify(ex, polymec_warn);
  assert_false(result);
  result = blob_exchanger_verify(ex, NULL);
  assert_false(result);

  // Now try to exchange data.
  int small_blob = 0, big_blob = 1;
  blob_buffer_t* buffer = blob_exchanger_create_buffer(ex, 1);
  small_widget_t smallw;
  blob_exchanger_copy_in(ex, small_blob, 1, &smallw, buffer);
  big_widget_t bigw;
  blob_exchanger_copy_in(ex, big_blob, 1, &bigw, buffer);
  blob_exchanger_exchange(ex, 0, buffer);

  // Put everything away.
  blob_buffer_free(buffer);
  release_ref(ex);

  // Next, we put together a good exchanger and check that verify/deadlock
  // detection give it a pass.
  ex = ring_exchanger(state);
  blob_exchanger_enable_deadlock_detection(ex, grace_period, 0, stdout);
  assert_true(blob_exchanger_deadlock_detection_enabled(ex));

  // Make sure the exchange works.
  buffer = blob_exchanger_create_buffer(ex, 1);
  blob_exchanger_copy_in(ex, small_blob, 1, &smallw, buffer);
  blob_exchanger_copy_in(ex, big_blob, 1, &bigw, buffer);
  blob_exchanger_exchange(ex, 0, buffer);

  // Clean up.
  blob_buffer_free(buffer);
  release_ref(ex);
}
*/

int main(int argc, char* argv[])
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] =
  {
    cmocka_unit_test(test_blob_exchanger_construct),
    cmocka_unit_test(test_blob_exchanger_exchange),
//    cmocka_unit_test(test_blob_exchanger_verify_and_dl_detection),
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
