// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "geometry/unimesh_patch.h"
#include "geometry/unimesh_patch_bc.h"

typedef struct
{
  real_t* values;
  int num_components;
} constant_bc_t;

static void update_x1_cells(void* context, unimesh_t* mesh,
                            int i, int j, int k, real_t t,
                            unimesh_patch_t* patch)
{
  constant_bc_t* bc = context;
  ASSERT(bc->num_components == patch->nc);

  DECLARE_UNIMESH_CELL_ARRAY(a, patch);
  for (int jj = 1; jj <= patch->ny; ++jj)
    for (int kk = 1; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[0][jj][kk][c] = bc->values[c];
}

static void update_x2_cells(void* context, unimesh_t* mesh,
                            int i, int j, int k, real_t t,
                            unimesh_patch_t* patch)
{
  constant_bc_t* bc = context;
  ASSERT(bc->num_components == patch->nc);

  DECLARE_UNIMESH_CELL_ARRAY(a, patch);
  for (int jj = 1; jj <= patch->ny; ++jj)
    for (int kk = 1; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[patch->nx+1][jj][kk][c] = bc->values[c];
}

static void update_y1_cells(void* context, unimesh_t* mesh,
                            int i, int j, int k, real_t t,
                            unimesh_patch_t* patch)
{
  constant_bc_t* bc = context;
  ASSERT(bc->num_components == patch->nc);

  DECLARE_UNIMESH_CELL_ARRAY(a, patch);
  for (int ii = 1; ii <= patch->nx; ++ii)
    for (int kk = 1; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][0][kk][c] = bc->values[c];
}

static void update_y2_cells(void* context, unimesh_t* mesh,
                            int i, int j, int k, real_t t,
                            unimesh_patch_t* patch)
{
  constant_bc_t* bc = context;
  ASSERT(bc->num_components == patch->nc);

  DECLARE_UNIMESH_CELL_ARRAY(a, patch);
  for (int ii = 1; ii <= patch->nx; ++ii)
    for (int kk = 1; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][patch->ny+1][kk][c] = bc->values[c];
}

static void update_z1_cells(void* context, unimesh_t* mesh,
                            int i, int j, int k, real_t t,
                            unimesh_patch_t* patch)
{
  constant_bc_t* bc = context;
  ASSERT(bc->num_components == patch->nc);

  DECLARE_UNIMESH_CELL_ARRAY(a, patch);
  for (int ii = 1; ii <= patch->nx; ++ii)
    for (int jj = 1; jj <= patch->ny; ++jj)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][jj][0][c] = bc->values[c];
}

static void update_z2_cells(void* context, unimesh_t* mesh,
                            int i, int j, int k, real_t t,
                            unimesh_patch_t* patch)
{
  constant_bc_t* bc = context;
  ASSERT(bc->num_components == patch->nc);

  DECLARE_UNIMESH_CELL_ARRAY(a, patch);
  for (int ii = 1; ii <= patch->nx; ++ii)
    for (int jj = 1; jj <= patch->ny; ++jj)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][jj][patch->nz+1][c] = bc->values[c];
}

static void update_x1_nodes(void* context, unimesh_t* mesh,
                            int i, int j, int k, real_t t,
                            unimesh_patch_t* patch)
{
  constant_bc_t* bc = context;
  ASSERT(bc->num_components == patch->nc);

  DECLARE_UNIMESH_NODE_ARRAY(a, patch);
  for (int jj = 0; jj <= patch->ny; ++jj)
    for (int kk = 0; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[0][jj][kk][c] = bc->values[c];
}

static void update_x2_nodes(void* context, unimesh_t* mesh,
                            int i, int j, int k, real_t t,
                            unimesh_patch_t* patch)
{
  constant_bc_t* bc = context;
  ASSERT(bc->num_components == patch->nc);

  DECLARE_UNIMESH_NODE_ARRAY(a, patch);
  for (int jj = 0; jj <= patch->ny; ++jj)
    for (int kk = 0; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[patch->nx][jj][kk][c] = bc->values[c];
}

static void update_y1_nodes(void* context, unimesh_t* mesh,
                            int i, int j, int k, real_t t,
                            unimesh_patch_t* patch)
{
  constant_bc_t* bc = context;
  ASSERT(bc->num_components == patch->nc);

  DECLARE_UNIMESH_NODE_ARRAY(a, patch);
  for (int ii = 0; ii <= patch->nx; ++ii)
    for (int kk = 0; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][0][kk][c] = bc->values[c];
}

static void update_y2_nodes(void* context, unimesh_t* mesh,
                            int i, int j, int k, real_t t,
                            unimesh_patch_t* patch)
{
  constant_bc_t* bc = context;
  ASSERT(bc->num_components == patch->nc);

  DECLARE_UNIMESH_NODE_ARRAY(a, patch);
  for (int ii = 0; ii <= patch->nx; ++ii)
    for (int kk = 0; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][patch->ny][kk][c] = bc->values[c];
}

static void update_z1_nodes(void* context, unimesh_t* mesh,
                            int i, int j, int k, real_t t,
                            unimesh_patch_t* patch)
{
  constant_bc_t* bc = context;
  ASSERT(bc->num_components == patch->nc);

  DECLARE_UNIMESH_NODE_ARRAY(a, patch);
  for (int ii = 0; ii <= patch->nx; ++ii)
    for (int jj = 0; jj <= patch->ny; ++jj)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][jj][0][c] = bc->values[c];
}

static void update_z2_nodes(void* context, unimesh_t* mesh,
                            int i, int j, int k, real_t t,
                            unimesh_patch_t* patch)
{
  constant_bc_t* bc = context;
  ASSERT(bc->num_components == patch->nc);

  DECLARE_UNIMESH_NODE_ARRAY(a, patch);
  for (int ii = 0; ii <= patch->nx; ++ii)
    for (int jj = 0; jj <= patch->ny; ++jj)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][jj][patch->nz][c] = bc->values[c];
}

static void constant_bc_free(void* context)
{
  constant_bc_t* bc = context;
  polymec_free(bc->values);
  polymec_free(bc);
}

unimesh_patch_bc_t* constant_unimesh_patch_bc_new(unimesh_t* mesh, 
                                                  real_t* values, 
                                                  int num_components)
{
  ASSERT(values != NULL);
  ASSERT(num_components > 0);

  constant_bc_t* bc = polymec_malloc(sizeof(constant_bc_t));
  bc->values = polymec_malloc(sizeof(real_t) * num_components);
  memcpy(bc->values, values, sizeof(real_t) * num_components);
  bc->num_components = num_components;

  unimesh_patch_bc_vtable vtable = {.dtor = constant_bc_free};
  vtable.start_update[0][0] = update_x1_cells;
  vtable.start_update[0][1] = update_x2_cells;
  vtable.start_update[0][2] = update_y1_cells;
  vtable.start_update[0][3] = update_y2_cells;
  vtable.start_update[0][4] = update_z1_cells;
  vtable.start_update[0][5] = update_z2_cells;
#if 0
  vtable.start_update[1][0] = start_update_xface_x1;
  vtable.start_update[1][1] = start_update_xface_x2;
  vtable.start_update[1][2] = start_update_xface_y1;
  vtable.start_update[1][3] = start_update_xface_y2;
  vtable.start_update[1][4] = start_update_xface_z1;
  vtable.start_update[1][5] = start_update_xface_z2;
  vtable.start_update[2][0] = start_update_yface_x1;
  vtable.start_update[2][1] = start_update_yface_x2;
  vtable.start_update[2][2] = start_update_yface_y1;
  vtable.start_update[2][3] = start_update_yface_y2;
  vtable.start_update[2][4] = start_update_yface_z1;
  vtable.start_update[2][5] = start_update_yface_z2;
  vtable.start_update[3][0] = start_update_zface_x1;
  vtable.start_update[3][1] = start_update_zface_x2;
  vtable.start_update[3][2] = start_update_zface_y1;
  vtable.start_update[3][3] = start_update_zface_y2;
  vtable.start_update[3][4] = start_update_zface_z1;
  vtable.start_update[3][5] = start_update_zface_z2;
  vtable.start_update[4][0] = start_update_xedge_x1;
  vtable.start_update[4][1] = start_update_xedge_x2;
  vtable.start_update[4][2] = start_update_xedge_y1;
  vtable.start_update[4][3] = start_update_xedge_y2;
  vtable.start_update[4][4] = start_update_xedge_z1;
  vtable.start_update[4][5] = start_update_xedge_z2;
  vtable.start_update[5][0] = start_update_yedge_x1;
  vtable.start_update[5][1] = start_update_yedge_x2;
  vtable.start_update[5][2] = start_update_yedge_y1;
  vtable.start_update[5][3] = start_update_yedge_y2;
  vtable.start_update[5][4] = start_update_yedge_z1;
  vtable.start_update[5][5] = start_update_yedge_z2;
  vtable.start_update[6][0] = start_update_zedge_x1;
  vtable.start_update[6][1] = start_update_zedge_x2;
  vtable.start_update[6][2] = start_update_zedge_y1;
  vtable.start_update[6][3] = start_update_zedge_y2;
  vtable.start_update[6][4] = start_update_zedge_z1;
  vtable.start_update[6][5] = start_update_zedge_z2;
#endif
  vtable.start_update[7][0] = update_x1_nodes;
  vtable.start_update[7][1] = update_x2_nodes;
  vtable.start_update[7][2] = update_y1_nodes;
  vtable.start_update[7][3] = update_y2_nodes;
  vtable.start_update[7][4] = update_z1_nodes;
  vtable.start_update[7][5] = update_z2_nodes;

  char name[257];
  snprintf(name, 256, "constant bc (%d components)", num_components);
  return unimesh_patch_bc_new(name, bc, vtable, mesh);
}
