// Copyright (c) 2012-2014, Jeffrey N. Johnson
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this 
// list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice, 
// this list of conditions and the following disclaimer in the documentation 
// and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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
  log_debug("preconditioner: setting up...");
  precond->vtable.setup(precond->context);
}

bool preconditioner_solve(preconditioner_t* precond, real_t* rhs)
{
  log_debug("preconditioner: solving...");
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

