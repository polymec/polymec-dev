// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "solvers/krylov_solver.h"

#if POLYMEC_HAVE_SHARED_LIBS
typedef int32_t HYPRE_Int;
#define HypreFactory hypre_krylov_factory_impl_32
#include "solvers/hypre_krylov_solver_impl.c"
#endif

krylov_factory_t* hypre_krylov_factory_32(const char* hypre_dir);
krylov_factory_t* hypre_krylov_factory_32(const char* hypre_dir)
{
#if POLYMEC_HAVE_SHARED_LIBS
  return HypreFactory(hypre_dir);
#else
  log_urgent("hypre_krylov_factory: Polymec must be configured with shared library support to use HYPRE.");
  return NULL;
#endif
}

