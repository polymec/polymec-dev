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
#include "simple_model.h"

static void test_ctor(void** state)
{
  model_t* m = simple_new();
  assert_true(strcmp(model_name(m), "simple") == 0);
  assert_true(model_interpreter(m) != NULL);
  assert_true(model_parallelism(m) == MODEL_MPI);
  assert_true(reals_equal(model_time(m), 0.0));
  assert_int_equal(0, model_step(m));
  model_free(m);
}

static void test_run(void** state)
{
  model_t* m = simple_new();
  model_read_input_string(m, "index = 0\ndt = 0.1\n");
  model_run(m, 0.0, 1.0, 10);
  assert_true((reals_nearly_equal(model_time(m), 1.0, 1e-12) || 
              (model_step(m) == 10)));
  model_free(m);

  // Check for the existence of the file.
  int nproc, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0)
  {
    char filename[FILENAME_MAX];
    snprintf(filename, FILENAME_MAX-1, "simple.%d.%d.%d", 0, nproc, nproc);
    assert_true(file_exists(filename));

    // Remove it.
    remove(filename);
  }
}

static void test_run_files(void** state)
{
  int rank, nproc;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  int num_files = 2*nproc;

  // Rank 0 writes input files with different time steps.
  char* inputs[num_files];
  for (int i = 0; i < num_files; ++i)
  {
    char filename[FILENAME_MAX];
    snprintf(filename, FILENAME_MAX-1, "simple_%d", i);
    if (rank == 0)
    {
      FILE* file = fopen(filename, "w");
      fprintf(file, "index = %d\n", i);
      fprintf(file, "dt = %g\n", (i+1) * 0.1);
      fprintf(file, "t2 = %g\n", 1.0 * (i+1));
      fclose(file);
    }
    inputs[i] = string_dup(filename);
  }

  // We will run each simulation with 1 or 2 processes.
  int procs_per_run = MIN(nproc, 2);

  // Now we create a model and run it on all the files.
  MPI_Barrier(MPI_COMM_WORLD);
  model_t* m = simple_new();
  model_run_files(m, inputs, num_files, procs_per_run);
  model_free(m);
  MPI_Barrier(MPI_COMM_WORLD);

  // Rank 0 checks for output and deletes the input/output files.
  if (rank == 0)
  {
    for (int i = 0; i < num_files; ++i)
    {
      // Input files.
      remove(inputs[i]);
      string_free(inputs[i]);

      // Output files.
      char filename[FILENAME_MAX];
      snprintf(filename, FILENAME_MAX-1, "simple.%d.%d.%d", i, nproc, procs_per_run);
      log_debug("Searching for output: %s", filename);
      assert_true(file_exists(filename));
      remove(filename);
    }
  }
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_ctor),
    cmocka_unit_test(test_run),
    cmocka_unit_test(test_run_files)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
