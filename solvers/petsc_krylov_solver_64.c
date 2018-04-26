// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "solvers/krylov_solver.h"

#if POLYMEC_HAVE_SHARED_LIBS
typedef int64_t PetscInt;
typedef real_t PetscScalar;
#define PetscFactory petsc_krylov_factory_impl_64
#include "solvers/petsc_krylov_solver_impl.c"
#endif

krylov_factory_t* petsc_krylov_factory_64(const char* petsc_library);
krylov_factory_t* petsc_krylov_factory_64(const char* petsc_library)
{
#if POLYMEC_HAVE_SHARED_LIBS
  return PetscFactory(petsc_library);
#else
  log_urgent("petsc_krylov_factory: Polymec must be configured with shared library support to use PETSc.");
  return NULL;
#endif
}

