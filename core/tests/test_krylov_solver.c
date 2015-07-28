// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>

#include "cmockery.h"
#include "core/polymec.h"
#include "core/krylov_solver.h"

void test_petsc_krylov_factory(void** state)
{
  char* petsc_dir = getenv("PETSC_DIR");
  if (petsc_dir != NULL)
  {
    char* petsc_arch = getenv("PETSC_ARCH");
    if (petsc_arch != NULL)
    {
      krylov_factory_t* petsc = petsc_krylov_factory(petsc_dir, petsc_arch);
      assert_true(petsc != NULL);
      assert_true(krylov_factory_name(petsc) != NULL);
      log_urgent("Created krylov_factory '%s'.", krylov_factory_name(petsc));
      krylov_factory_free(petsc);
    }
    else
      log_urgent("PETSC_ARCH not set. Skipping test_petsc_krylov_factory.");
  }
  else
    log_urgent("PETSC_DIR not set. Skipping test_petsc_krylov_factory.");
}

void test_petsc_krylov_vector(void** state)
{
  char* petsc_dir = getenv("PETSC_DIR");
  if (petsc_dir != NULL)
  {
    char* petsc_arch = getenv("PETSC_ARCH");
    if (petsc_arch != NULL)
    {
      krylov_factory_t* petsc = petsc_krylov_factory(petsc_dir, petsc_arch);
      krylov_vector_t* vec = krylov_factory_vector(petsc, MPI_COMM_WORLD, 100);
      krylov_vector_free(vec);
      krylov_factory_free(petsc);
    }
    else
      log_urgent("PETSC_ARCH not set. Skipping test_petsc_krylov_vector.");
  }
  else
    log_urgent("PETSC_DIR not set. Skipping test_petsc_krylov_vector.");
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_petsc_krylov_factory),
    unit_test(test_petsc_krylov_vector)
  };
  return run_tests(tests);
}

