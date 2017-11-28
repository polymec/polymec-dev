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
#include "core/unordered_set.h"
#include "core/timer.h"
#include "geometry/cubic_lattice.h"
#include "geometry/create_uniform_polymesh.h"
#include "model/polymesh_stencils.h"
#include "io/silo_file.h"

static void test_NXxNYxNZ_star_stencil(void** state, 
                                       MPI_Comm comm,
                                       int radius,
                                       int nx, int ny, int nz,
                                       int num_interior_neighbors, 
                                       int num_boundary_neighbors, 
                                       int num_edge_neighbors,
                                       int num_corner_neighbors)
{
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  polymesh_t* mesh = create_uniform_polymesh(comm, nx, ny, nz, &bbox);
  stencil_t* stencil = cell_star_stencil_new(mesh, radius);

  // Write the stencil to a Silo file.
  char prefix[FILENAME_MAX+1];
  if (comm == MPI_COMM_WORLD)
    snprintf(prefix, FILENAME_MAX, "star_all");
  else
  {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    snprintf(prefix, FILENAME_MAX, "star_p%d", rank);
  }
  silo_file_t* silo = silo_file_new(comm, prefix, "star", 1, 0, 0, 0.0);
  silo_file_write_stencil(silo, "stencil", stencil);
  silo_file_close(silo);

  // Make sure the file is written before we attempt to read it!
  MPI_Barrier(comm);

  // Read the stencil from the file and check its contents.
  real_t t;
  silo = silo_file_open(comm, prefix, "star", 0, 0, &t);
  assert_true(silo_file_contains_stencil(silo, "stencil"));
  stencil_t* stencil1 = silo_file_read_stencil(silo, "stencil", comm);
  silo_file_close(silo);
  int N = stencil_num_indices(stencil);
  assert_int_equal(N, stencil_num_indices(stencil1));
  for (int i = 0; i < N; ++i)
  {
    int pos = 0, pos1 = 0, j, j1;
    while (stencil_next(stencil, i, &pos, &j))
    {
      assert_true(stencil_next(stencil1, i, &pos1, &j1));
      assert_int_equal(j, j1);
    }
  }
  stencil_free(stencil1);
  stencil_free(stencil);
  polymesh_free(mesh);
}

static void test_serial_1x1x1_cell_star_stencil(void** state)
{
  START_FUNCTION_TIMER();
  test_NXxNYxNZ_star_stencil(state, MPI_COMM_SELF, 1, 1, 1, 1, 0, 0, 0, 0);
  STOP_FUNCTION_TIMER();
}

static void test_serial_10x1x1_cell_star_stencil(void** state)
{
  START_FUNCTION_TIMER();
  test_NXxNYxNZ_star_stencil(state, MPI_COMM_SELF, 1, 10, 1, 1, 0, 0, 2, 1);
  STOP_FUNCTION_TIMER();
}

static void test_serial_10x10x1_cell_star_stencil(void** state)
{
  START_FUNCTION_TIMER();
  test_NXxNYxNZ_star_stencil(state, MPI_COMM_SELF, 1, 10, 10, 1, 4, 4, 3, 2);
  STOP_FUNCTION_TIMER();
}

static void test_serial_10x10x10_cell_star_stencil(void** state)
{
  test_NXxNYxNZ_star_stencil(state, MPI_COMM_SELF, 1, 10, 10, 10, 6, 5, 4, 3);
  START_FUNCTION_TIMER();
  STOP_FUNCTION_TIMER();
}

#if POLYMEC_HAVE_MPI
static void test_parallel_10x1x1_cell_star_stencil(void** state)
{
  START_FUNCTION_TIMER();
  test_NXxNYxNZ_star_stencil(state, MPI_COMM_WORLD, 1, 10, 1, 1, 0, 0, 2, 1);
  STOP_FUNCTION_TIMER();
}

static void test_parallel_10x10x1_cell_star_stencil(void** state)
{
  START_FUNCTION_TIMER();
  test_NXxNYxNZ_star_stencil(state, MPI_COMM_WORLD, 1, 10, 10, 1, 4, 4, 3, 2);
  STOP_FUNCTION_TIMER();
}

static void test_parallel_10x10x10_cell_star_stencil(void** state)
{
  START_FUNCTION_TIMER();
  test_NXxNYxNZ_star_stencil(state, MPI_COMM_WORLD, 1, 10, 10, 10, 6, 5, 4, 3);
  STOP_FUNCTION_TIMER();
}
#endif

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_serial_1x1x1_cell_star_stencil),
    cmocka_unit_test(test_serial_10x1x1_cell_star_stencil),
    cmocka_unit_test(test_serial_10x10x1_cell_star_stencil),
    cmocka_unit_test(test_serial_10x10x10_cell_star_stencil)
#if POLYMEC_HAVE_MPI
   ,cmocka_unit_test(test_parallel_10x1x1_cell_star_stencil),
    cmocka_unit_test(test_parallel_10x10x1_cell_star_stencil),
    cmocka_unit_test(test_parallel_10x10x10_cell_star_stencil)
#endif
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
