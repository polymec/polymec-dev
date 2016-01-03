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
#include "core/file_utils.h"
#include "core/krylov_solver.h"

static krylov_factory_t* create_petsc_krylov_factory()
{
  krylov_factory_t* factory = NULL;
  char* petsc_dir = getenv("PETSC_DIR");
  if (petsc_dir != NULL)
  {
    char* petsc_arch = getenv("PETSC_ARCH");
    if (petsc_arch != NULL)
    {
      // Check for the existence of the PETSc dynamic library.
      char petsc_path[FILENAME_MAX+1];
      snprintf(petsc_path, FILENAME_MAX, "%s/%s/lib/libpetsc%s", petsc_dir, petsc_arch, SHARED_LIBRARY_SUFFIX);
      if (file_exists(petsc_path))
      {
        factory = petsc_krylov_factory(petsc_dir, petsc_arch);
        if (factory == NULL)
          log_urgent("Could not load PETSc. Skipping PETSc test.");
      }
      else
        log_urgent("PETSc library not found. Skipping PETSc test.");
    }
    else
      log_urgent("PETSC_ARCH not set. Skipping PETSc test.");
  }
  else
    log_urgent("PETSC_DIR not set. Skipping PETSc test.");
  return factory;
}

void test_petsc_krylov_factory(void** state)
{
  krylov_factory_t* petsc = create_petsc_krylov_factory();
  if (petsc != NULL)
  {
    assert_true(krylov_factory_name(petsc) != NULL);
    krylov_factory_free(petsc);
  }
}

void test_petsc_krylov_vector(void** state)
{
  krylov_factory_t* petsc = create_petsc_krylov_factory();
  if (petsc != NULL)
  {
    // Create a vector.
    krylov_vector_t* vec = krylov_factory_vector(petsc, MPI_COMM_WORLD, 100);
    assert_int_equal(100, krylov_vector_local_size(vec));
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    assert_int_equal(100*nprocs, krylov_vector_global_size(vec));

    // Clone it.
    krylov_vector_t* vec1 = krylov_vector_clone(vec);
    assert_int_equal(100, krylov_vector_local_size(vec1));

    // Put everything away.
    krylov_vector_free(vec);
    krylov_vector_free(vec1);
    krylov_factory_free(petsc);
  }
}

static krylov_factory_t* create_hypre_krylov_factory()
{
  krylov_factory_t* factory = NULL;
  char* hypre_dir = getenv("HYPRE_DIR"); // Absent other info, we rely on this.
  if (hypre_dir != NULL)
  {
    // Check for the existence of the HYPRE dynamic library.
    char hypre_path[FILENAME_MAX+1];
    snprintf(hypre_path, FILENAME_MAX, "%s/libHYPRE%s", hypre_dir, SHARED_LIBRARY_SUFFIX);
    if (file_exists(hypre_path))
    {
      factory = hypre_krylov_factory(hypre_dir);
      if (factory == NULL)
        log_urgent("Could not load HYPRE. Skipping HYPRE test.");
    }
    else
      log_urgent("HYPRE library not found. Skipping HYPRE test.");
  }
  else
    log_urgent("HYPRE_DIR not set. Skipping HYPRE test.");
  return factory;
}

void test_hypre_krylov_factory(void** state)
{
  krylov_factory_t* hypre = create_hypre_krylov_factory();
  if (hypre != NULL)
  {
    assert_true(krylov_factory_name(hypre) != NULL);
    krylov_factory_free(hypre);
  }
}


int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_petsc_krylov_factory),
    unit_test(test_petsc_krylov_vector),
    unit_test(test_hypre_krylov_factory)
  };
  return run_tests(tests);
}

