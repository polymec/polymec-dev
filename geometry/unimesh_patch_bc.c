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

unimesh_t* unimesh_patch_bc_mesh(unimesh_patch_bc_t* bc)
{
  return bc->mesh;
}

bool unimesh_patch_bc_handles_centering(unimesh_patch_bc_t* bc,
                                        unimesh_centering_t centering)
{
  return (((centering == UNIMESH_CELL) && (bc->vtable.start_cell_update != NULL)) ||
          ((centering == UNIMESH_XFACE) && (bc->vtable.start_xface_update != NULL)) ||
          ((centering == UNIMESH_YFACE) && (bc->vtable.start_yface_update != NULL)) ||
          ((centering == UNIMESH_ZFACE) && (bc->vtable.start_zface_update != NULL)) ||
          ((centering == UNIMESH_XEDGE) && (bc->vtable.start_xedge_update != NULL)) ||
          ((centering == UNIMESH_YEDGE) && (bc->vtable.start_yedge_update != NULL)) ||
          ((centering == UNIMESH_ZEDGE) && (bc->vtable.start_zedge_update != NULL)) ||
          ((centering == UNIMESH_NODE) && (bc->vtable.start_node_update != NULL)));
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
  void (*start_update)(void*, unimesh_t*, int, int, int, real_t, 
                       unimesh_boundary_t, unimesh_patch_t*);
  switch(patch->centering) 
  {
    case UNIMESH_CELL:
      start_update = bc->vtable.start_cell_update; break;
    case UNIMESH_XFACE:
      start_update = bc->vtable.start_xface_update; break;
    case UNIMESH_YFACE:
      start_update = bc->vtable.start_yface_update; break;
    case UNIMESH_ZFACE:
      start_update = bc->vtable.start_zface_update; break;
    case UNIMESH_XEDGE:
      start_update = bc->vtable.start_xedge_update; break;
    case UNIMESH_YEDGE:
      start_update = bc->vtable.start_yedge_update; break;
    case UNIMESH_ZEDGE:
      start_update = bc->vtable.start_zedge_update; break;
    default:
      start_update = bc->vtable.start_node_update; 
  }
  ASSERT(start_update != NULL);
  start_update(bc->context, bc->mesh, i, j, k, t, patch_boundary, patch);
}

void unimesh_patch_bc_finish_update(unimesh_patch_bc_t* bc,
                                    int i, int j, int k, real_t t,
                                    unimesh_boundary_t patch_boundary,
                                    unimesh_patch_t* patch)
{
  void (*finish_update)(void*, unimesh_t*, int, int, int, real_t, 
                        unimesh_boundary_t, unimesh_patch_t*);
  switch(patch->centering) 
  {
    case UNIMESH_CELL:
      finish_update = bc->vtable.finish_cell_update; break;
    case UNIMESH_XFACE:
      finish_update = bc->vtable.finish_xface_update; break;
    case UNIMESH_YFACE:
      finish_update = bc->vtable.finish_yface_update; break;
    case UNIMESH_ZFACE:
      finish_update = bc->vtable.finish_zface_update; break;
    case UNIMESH_XEDGE:
      finish_update = bc->vtable.finish_xedge_update; break;
    case UNIMESH_YEDGE:
      finish_update = bc->vtable.finish_yedge_update; break;
    case UNIMESH_ZEDGE:
      finish_update = bc->vtable.finish_zedge_update; break;
    default:
      finish_update = bc->vtable.finish_node_update;
  }
  if (finish_update != NULL)
    finish_update(bc->context, bc->mesh, i, j, k, t, patch_boundary, patch);
}

