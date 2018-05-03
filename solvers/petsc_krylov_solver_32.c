// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "solvers/krylov_solver.h"

typedef int32_t PetscInt;
typedef real_t PetscScalar;
#define PetscFactory petsc_krylov_factory_impl_32
#include "solvers/petsc_krylov_solver_impl.c"

krylov_factory_t* petsc_krylov_factory_32(const char* petsc_library);
krylov_factory_t* petsc_krylov_factory_32(const char* petsc_library)
{
  return PetscFactory(petsc_library);
}

