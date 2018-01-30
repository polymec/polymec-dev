// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "geometry/unimesh_patch_bc.h"

struct unimesh_patch_bc_t 
{
  char* name;
  void* context;
  unimesh_patch_bc_vtable vtable;
  unimesh_t* mesh;

  int num_components;
  int* components;
};

static void unimesh_patch_bc_free(void* context)
{
  unimesh_patch_bc_t* bc = context;
  if ((bc->context != NULL) && (bc->vtable.dtor != NULL))
    bc->vtable.dtor(bc->context);
  string_free(bc->name);
  polymec_free(bc->components);
}

unimesh_patch_bc_t* unimesh_patch_bc_new(const char* name,
                                         void* context,
                                         unimesh_patch_bc_vtable vtable,
                                         unimesh_t* mesh)
{
  ASSERT(vtable.handles_centering != NULL);
  ASSERT(vtable.start_update != NULL);
  ASSERT(vtable.finish_update != NULL);

  unimesh_patch_bc_t* bc = polymec_gc_malloc(sizeof(unimesh_patch_bc_t), 
                                             unimesh_patch_bc_free);
  bc->name = string_dup(name);
  bc->context = context;
  bc->vtable = vtable;
  bc->mesh = mesh;
  bc->num_components = 0; // all components
  return bc;
}

unimesh_patch_bc_t* multicomp_unimesh_patch_bc_new(const char* name,
                                                   void* context,
                                                   unimesh_patch_bc_vtable vtable,
                                                   unimesh_t* mesh,
                                                   int num_components)
{
  ASSERT(vtable.handles_centering != NULL);
  ASSERT(vtable.start_update != NULL);
  ASSERT(vtable.finish_update != NULL);
  ASSERT(num_components > 0);

  unimesh_patch_bc_t* bc = polymec_gc_malloc(sizeof(unimesh_patch_bc_t), 
                                             unimesh_patch_bc_free);
  bc->name = string_dup(name);
  bc->context = context;
  bc->vtable = vtable;
  bc->mesh = mesh;
  bc->num_components = num_components;
  return bc;
}

char* unimesh_patch_bc_name(unimesh_patch_bc_t* bc)
{
  return bc->name;
}

int unimesh_patch_bc_num_components(unimesh_patch_bc_t* bc)
{
  return bc->num_components;
}

unimesh_t* unimesh_patch_bc_mesh(unimesh_patch_bc_t* bc)
{
  return bc->mesh;
}

bool unimesh_patch_bc_handles_centering(unimesh_patch_bc_t* bc,
                                        unimesh_centering_t centering)
{
  return bc->vtable.handles_centering(bc->context, centering);
}

void unimesh_patch_bc_update(unimesh_patch_bc_t* bc,
                             int i, int j, int k, real_t t,
                             unimesh_boundary_t patch_boundary,
                             unimesh_patch_t* patch)
{
  unimesh_patch_bc_start_update(bc, i, j, k, t, patch_boundary, patch);
  unimesh_patch_bc_finish_update(bc, i, j, k, t, patch_boundary, patch);
}

void unimesh_patch_bc_start_update(unimesh_patch_bc_t* bc,
                                   int i, int j, int k, real_t t,
                                   unimesh_boundary_t patch_boundary,
                                   unimesh_patch_t* patch)
{
  bc->vtable.start_update(bc->context, bc->mesh, i, j, k, t,
                          bc->components, bc->num_components, 
                          patch_boundary, patch);
}

void unimesh_patch_bc_finish_update(unimesh_patch_bc_t* bc,
                                    int i, int j, int k, real_t t,
                                    unimesh_boundary_t patch_boundary,
                                    unimesh_patch_t* patch)
{
  bc->vtable.finish_update(bc->context, bc->mesh, i, j, k, t,
                           bc->components, bc->num_components, 
                           patch_boundary, patch);
}

