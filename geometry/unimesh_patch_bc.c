// Copyright (c) 2012-2019, Jeffrey N. Johnson
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
  unimesh_t* mesh;
  void* context;
  unimesh_patch_bc_vtable vtable;

  // Easy version stuff.
  void* easy_context;
  unimesh_patch_bc_easy_vtable easy_vtable;
};

static void unimesh_patch_bc_free(void* context)
{
  unimesh_patch_bc_t* bc = context;
  if ((bc->context != NULL) && (bc->vtable.dtor != NULL))
    bc->vtable.dtor(bc->context);
  string_free(bc->name);
}

unimesh_patch_bc_t* unimesh_patch_bc_new(const char* name,
                                         void* context,
                                         unimesh_patch_bc_vtable vtable,
                                         unimesh_t* mesh)
{
  unimesh_patch_bc_t* bc = polymec_refcounted_malloc(sizeof(unimesh_patch_bc_t),
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

  bc->easy_context = NULL;
  memset(&(bc->easy_vtable), 0, sizeof(unimesh_patch_bc_easy_vtable));
  bc->mesh = mesh;
  return bc;
}

char* unimesh_patch_bc_name(unimesh_patch_bc_t* bc)
{
  return bc->name;
}

void* unimesh_patch_bc_context(unimesh_patch_bc_t* bc)
{
  return (bc != NULL) ? (bc->easy_context != NULL) ? bc->easy_context : bc->context
                      : NULL;
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
                             field_metadata_t* md,
                             unimesh_patch_t* patch)
{
  unimesh_patch_bc_start_update(bc, i, j, k, t, patch_boundary, md, patch);
  unimesh_patch_bc_finish_update(bc, i, j, k, t, patch_boundary, md, patch);
}

void unimesh_patch_bc_start_update(unimesh_patch_bc_t* bc,
                                   int i, int j, int k, real_t t,
                                   unimesh_boundary_t patch_boundary,
                                   field_metadata_t* md,
                                   unimesh_patch_t* patch)
{
  ASSERT(unimesh_patch_bc_handles_centering(bc, patch->centering));

  int cent = (int)patch->centering;
  int bnd = (int)patch_boundary;
  unimesh_patch_bc_update_method start_update = bc->vtable.start_update[cent][bnd];
  start_update(bc->context, bc->mesh, i, j, k, t, md, patch);
}

void unimesh_patch_bc_finish_update(unimesh_patch_bc_t* bc,
                                    int i, int j, int k, real_t t,
                                    unimesh_boundary_t patch_boundary,
                                    field_metadata_t* md,
                                    unimesh_patch_t* patch)
{
  ASSERT(unimesh_patch_bc_handles_centering(bc, patch->centering));

  int cent = (int)patch->centering;
  int bnd = (int)patch_boundary;
  unimesh_patch_bc_update_method finish_update = bc->vtable.finish_update[cent][bnd];
  if (finish_update != NULL)
    finish_update(bc->context, bc->mesh, i, j, k, t, md, patch);
}

#define DEFINE_EASY_UPDATES(boundary_name, boundary) \
static void easy_##boundary_name##_start_update(void* context, \
                                                unimesh_t* mesh, \
                                                int i, int j, int k, real_t t, \
                                                field_metadata_t* md, \
                                                unimesh_patch_t* patch) \
{ \
  unimesh_patch_bc_t* bc = context; \
  bc->easy_vtable.start_update(bc->easy_context, mesh, i, j, k, t, boundary, md, patch); \
} \
\
static void easy_##boundary_name##_finish_update(void* context, \
                                                 unimesh_t* mesh, \
                                                 int i, int j, int k, real_t t, \
                                                 field_metadata_t* md, \
                                                 unimesh_patch_t* patch) \
{ \
  unimesh_patch_bc_t* bc = context; \
  if (bc->easy_vtable.finish_update != NULL) \
    bc->easy_vtable.finish_update(bc->easy_context, mesh, i, j, k, t, boundary, md, patch); \
} \

DEFINE_EASY_UPDATES(x1, UNIMESH_X1_BOUNDARY)
DEFINE_EASY_UPDATES(x2, UNIMESH_X2_BOUNDARY)
DEFINE_EASY_UPDATES(y1, UNIMESH_Y1_BOUNDARY)
DEFINE_EASY_UPDATES(y2, UNIMESH_Y2_BOUNDARY)
DEFINE_EASY_UPDATES(z1, UNIMESH_Z1_BOUNDARY)
DEFINE_EASY_UPDATES(z2, UNIMESH_Z2_BOUNDARY)

static void easy_dtor(void* context)
{
  unimesh_patch_bc_t* bc = context;
  if ((bc->easy_context != NULL) && (bc->easy_vtable.dtor != NULL))
    bc->easy_vtable.dtor(bc->easy_context);
}

unimesh_patch_bc_t* unimesh_patch_bc_new_easy(const char* name,
                                              void* context,
                                              unimesh_patch_bc_easy_vtable vtable,
                                              unimesh_t* mesh)
{
  unimesh_patch_bc_t* bc = polymec_refcounted_malloc(sizeof(unimesh_patch_bc_t),
                                                     unimesh_patch_bc_free);
  bc->name = string_dup(name);
  bc->context = bc;
  for (int c = 0; c < 8; ++c)
  {
    bc->vtable.start_update[c][0] = easy_x1_start_update;
    bc->vtable.start_update[c][1] = easy_x2_start_update;
    bc->vtable.start_update[c][2] = easy_y1_start_update;
    bc->vtable.start_update[c][3] = easy_y2_start_update;
    bc->vtable.start_update[c][4] = easy_z1_start_update;
    bc->vtable.start_update[c][5] = easy_z2_start_update;

    bc->vtable.finish_update[c][0] = easy_x1_finish_update;
    bc->vtable.finish_update[c][1] = easy_x2_finish_update;
    bc->vtable.finish_update[c][2] = easy_y1_finish_update;
    bc->vtable.finish_update[c][3] = easy_y2_finish_update;
    bc->vtable.finish_update[c][4] = easy_z1_finish_update;
    bc->vtable.finish_update[c][5] = easy_z2_finish_update;
  }
  bc->vtable.dtor = easy_dtor;

  bc->easy_context = context;
  bc->easy_vtable = vtable;
  bc->mesh = mesh;
  return bc;
}

