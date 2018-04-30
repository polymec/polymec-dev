// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/options.h"
#include "core/string_utils.h"
#include "solvers/krylov_solver.h"

// These functions seamlessly create factories for PETSc- and HYPRE-based
// Krylov solvers.
krylov_factory_t* create_petsc_krylov_factory(void);
krylov_factory_t* create_petsc_krylov_factory()
{
  char petsc_library[FILENAME_MAX+1];
  petsc_library[0] = '\0';
  options_t* options = options_argv();
  char* petsc_lib = options_value(options, "petsc_library");
  char* petsc_64 = options_value(options, "petsc_64");
  bool use_64_bit_indices = false;
  if (petsc_64 != NULL)
    use_64_bit_indices = string_as_boolean(petsc_64);
  if (petsc_lib != NULL)
    strcpy(petsc_library, petsc_lib);
  else
  {
    log_urgent("petsc_library option not given. Using PETSC_DIR/PETSC_ARCH.");
    char* petsc_dir = getenv("PETSC_DIR");
    if (petsc_dir != NULL)
    {
      char* petsc_arch = getenv("PETSC_ARCH");
      if (petsc_arch != NULL)
        snprintf(petsc_library, FILENAME_MAX, "%s/%s/lib/libpetsc%s", petsc_dir, petsc_arch, SHARED_LIBRARY_SUFFIX);
      else
        log_urgent("PETSC_ARCH not set. Skipping PETSc test.");
    }
    else
      log_urgent("PETSC_DIR not set. Skipping PETSc test.");
  }
  krylov_factory_t* factory = NULL;
  if (strlen(petsc_library) > 0)
  {
    factory = petsc_krylov_factory(petsc_library, use_64_bit_indices);
    if (factory == NULL)
      log_urgent("Could not load PETSc. Skipping PETSc test.");
  }
  return factory;
}

krylov_factory_t* create_hypre_krylov_factory(void);
krylov_factory_t* create_hypre_krylov_factory()
{
  char hypre_dir[FILENAME_MAX+1];
  hypre_dir[0] = '\0';
  options_t* options = options_argv();
  char* hypre_opt_dir = options_value(options, "hypre_dir");
  char* hypre_64 = options_value(options, "hypre_64");
  bool use_64_bit_indices = false;
  if (hypre_64 != NULL)
    use_64_bit_indices = string_as_boolean(hypre_64);
  if (hypre_opt_dir != NULL)
    strcpy(hypre_dir, hypre_opt_dir);
  else
  {
    log_urgent("hypre_dir option not given. Using HYPRE_DIR.");
    char* hypre_env_dir = getenv("HYPRE_DIR"); // Absent other info, we rely on this.
    if (hypre_env_dir != NULL)
      strcpy(hypre_dir, hypre_env_dir);
    else
      log_urgent("HYPRE_DIR not set. Skipping HYPRE test.");
  }
  krylov_factory_t* factory = NULL;
  if (strlen(hypre_dir) > 0)
  {
    factory = hypre_krylov_factory(hypre_dir, use_64_bit_indices);
    if (factory == NULL)
      log_urgent("Could not load HYPRE. Skipping HYPRE test.");
  }
  return factory;
}


