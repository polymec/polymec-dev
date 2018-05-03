// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "solvers/krylov_solver.h"

typedef long long HYPRE_Int;
#define HypreFactory hypre_krylov_factory_impl_64
#include "solvers/hypre_krylov_solver_impl.c"

krylov_factory_t* hypre_krylov_factory_64(const char* hypre_dir);
krylov_factory_t* hypre_krylov_factory_64(const char* hypre_dir)
{
  return HypreFactory(hypre_dir);
}

