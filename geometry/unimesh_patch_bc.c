// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "geometry/unimesh_patch_bc.h"
#include "geometry/unimesh_patch.h"

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
  unimesh_patch_bc_t* bc = polymec_gc_malloc(sizeof(unimesh_patch_bc_t), 
                                             unimesh_patch_bc_free);
  bc->name = string_dup(name);
  bc->context = context;
  for (int c = 0; c < 8; ++c)
  {
    for (int b = 0; b < 6; ++b)
    {
      bc->vtable.start_update[c][b] = vtable.start_update[c][b];
      bc->vtable.finish_update[c][b] = vtable.finish_update[c][b];
    }
  }
  bc->vtable.dtor = vtable.dtor;
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
  ASSERT(num_components > 0);
  unimesh_patch_bc_t* bc = unimesh_patch_bc_new(name, context, vtable, mesh);
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

void* unimesh_patch_bc_context(unimesh_patch_bc_t* bc)
{
  return bc->context;
}

unimesh_t* unimesh_patch_bc_mesh(unimesh_patch_bc_t* bc)
{
  return bc->mesh;
}

bool unimesh_patch_bc_handles_centering(unimesh_patch_bc_t* bc,
                                        unimesh_centering_t centering)
{
  bool handles = true;
  int cent = (int)centering;
  for (int b = 0; b < 6; ++b)
  {
    if (bc->vtable.start_update[cent][b] == NULL)
    {
      handles = false; 
      break;
    }
  }
  return handles;
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
  ASSERT(unimesh_patch_bc_handles_centering(bc, patch->centering));

  int cent = (int)patch->centering;
  int bnd = (int)patch_boundary;
  unimesh_patch_bc_update_method start_update = bc->vtable.start_update[cent][bnd];
  start_update(bc->context, bc->mesh, i, j, k, t, patch);
}

void unimesh_patch_bc_finish_update(unimesh_patch_bc_t* bc,
                                    int i, int j, int k, real_t t,
                                    unimesh_boundary_t patch_boundary,
                                    unimesh_patch_t* patch)
{
  ASSERT(unimesh_patch_bc_handles_centering(bc, patch->centering));

  int cent = (int)patch->centering;
  int bnd = (int)patch_boundary;
  unimesh_patch_bc_update_method finish_update = bc->vtable.finish_update[cent][bnd];
  if (finish_update != NULL)
    finish_update(bc->context, bc->mesh, i, j, k, t, patch);
}

