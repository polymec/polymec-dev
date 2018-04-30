// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "solvers/krylov_solver.h"

extern krylov_factory_t* hypre_krylov_factory_32(const char* hypre_dir);
extern krylov_factory_t* hypre_krylov_factory_64(const char* hypre_dir);

krylov_factory_t* hypre_krylov_factory(const char* hypre_dir,
                                       bool use_64_bit_indices)
{
  if (!directory_exists(hypre_dir))
  {
    log_urgent("HYPRE library directory %s not found.", hypre_dir);
    return NULL;
  }
  if (use_64_bit_indices)
    return hypre_krylov_factory_64(hypre_dir);
  else
    return hypre_krylov_factory_32(hypre_dir);
}

