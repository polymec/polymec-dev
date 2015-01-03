// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/preconditioner.h"

struct preconditioner_t 
{
  char* name;
  void* context;
  preconditioner_vtable vtable;
};

preconditioner_t* preconditioner_new(const char* name,
                                     void* context,
                                     preconditioner_vtable vtable)
{
  ASSERT(vtable.setup != NULL);
  ASSERT(vtable.solve != NULL);
  preconditioner_t* precond = polymec_malloc(sizeof(preconditioner_t));
  precond->name = string_dup(name);
  precond->context = context;
  precond->vtable = vtable;

  return precond;
}

void preconditioner_free(preconditioner_t* precond)
{
  if ((precond->vtable.dtor != NULL) && (precond->context != NULL))
    precond->vtable.dtor(precond->context);
  polymec_free(precond->name);
  polymec_free(precond);
}

char* preconditioner_name(preconditioner_t* precond)
{
  return precond->name;
}

void* preconditioner_context(preconditioner_t* precond)
{
  return precond->context;
}

void preconditioner_setup(preconditioner_t* precond)
{
  precond->vtable.setup(precond->context);
}

bool preconditioner_solve(preconditioner_t* precond, real_t* rhs)
{
  return precond->vtable.solve(precond->context, rhs);
}

void preconditioner_fprintf(preconditioner_t* precond, FILE* stream)
{
  if (stream == NULL) return;
  if (precond->vtable.fprintf != NULL)
    precond->vtable.fprintf(precond->context, stream);
  else
    fprintf(stream, "Preconditioner '%s'", precond->name);
}

