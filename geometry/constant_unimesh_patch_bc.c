// Copyright (c) 2012-2019, Jeffrey N. Johnson
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
                            field_metadata_t* md,
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
                            field_metadata_t* md,
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
                            field_metadata_t* md,
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
                            field_metadata_t* md,
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
                            field_metadata_t* md,
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
                            field_metadata_t* md,
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

static void update_x1_xfaces(void* context, unimesh_t* mesh,
                             int i, int j, int k, real_t t,
                             field_metadata_t* md,
                             unimesh_patch_t* patch)
{
  constant_bc_t* bc = context;
  ASSERT(bc->num_components == patch->nc);

  DECLARE_UNIMESH_XFACE_ARRAY(a, patch);
  for (int jj = 0; jj < patch->ny; ++jj)
    for (int kk = 0; kk < patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[0][jj][kk][c] = bc->values[c];
}

static void update_x2_xfaces(void* context, unimesh_t* mesh,
                             int i, int j, int k, real_t t,
                             field_metadata_t* md,
                             unimesh_patch_t* patch)
{
  constant_bc_t* bc = context;
  ASSERT(bc->num_components == patch->nc);

  DECLARE_UNIMESH_XFACE_ARRAY(a, patch);
  for (int jj = 0; jj < patch->ny; ++jj)
    for (int kk = 0; kk < patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[patch->nx][jj][kk][c] = bc->values[c];
}

static void update_y1_xfaces(void* context, unimesh_t* mesh,
                             int i, int j, int k, real_t t,
                             field_metadata_t* md,
                             unimesh_patch_t* patch)
{
  // X faces aren't communicated across y boundaries.
}

static void update_y2_xfaces(void* context, unimesh_t* mesh,
                             int i, int j, int k, real_t t,
                             field_metadata_t* md,
                             unimesh_patch_t* patch)
{
  // X faces aren't communicated across y boundaries.
}

static void update_z1_xfaces(void* context, unimesh_t* mesh,
                             int i, int j, int k, real_t t,
                             field_metadata_t* md,
                             unimesh_patch_t* patch)
{
  // X faces aren't communicated across z boundaries.
}

static void update_z2_xfaces(void* context, unimesh_t* mesh,
                             int i, int j, int k, real_t t,
                             field_metadata_t* md,
                             unimesh_patch_t* patch)
{
  // X faces aren't communicated across z boundaries.
}

static void update_x1_yfaces(void* context, unimesh_t* mesh,
                             int i, int j, int k, real_t t,
                             field_metadata_t* md,
                             unimesh_patch_t* patch)
{
  // Y faces aren't communicated across x boundaries.
}

static void update_x2_yfaces(void* context, unimesh_t* mesh,
                             int i, int j, int k, real_t t,
                             field_metadata_t* md,
                             unimesh_patch_t* patch)
{
  // Y faces aren't communicated across x boundaries.
}

static void update_y1_yfaces(void* context, unimesh_t* mesh,
                             int i, int j, int k, real_t t,
                             field_metadata_t* md,
                             unimesh_patch_t* patch)
{
  constant_bc_t* bc = context;
  ASSERT(bc->num_components == patch->nc);

  DECLARE_UNIMESH_YFACE_ARRAY(a, patch);
  for (int ii = 0; ii < patch->nx; ++ii)
    for (int kk = 0; kk < patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][0][kk][c] = bc->values[c];
}

static void update_y2_yfaces(void* context, unimesh_t* mesh,
                             int i, int j, int k, real_t t,
                             field_metadata_t* md,
                             unimesh_patch_t* patch)
{
  constant_bc_t* bc = context;
  ASSERT(bc->num_components == patch->nc);

  DECLARE_UNIMESH_YFACE_ARRAY(a, patch);
  for (int ii = 0; ii < patch->nx; ++ii)
    for (int kk = 0; kk < patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][patch->ny][kk][c] = bc->values[c];
}

static void update_z1_yfaces(void* context, unimesh_t* mesh,
                             int i, int j, int k, real_t t,
                             field_metadata_t* md,
                             unimesh_patch_t* patch)
{
  // Y faces aren't communicated across z boundaries.
}

static void update_z2_yfaces(void* context, unimesh_t* mesh,
                             int i, int j, int k, real_t t,
                             field_metadata_t* md,
                             unimesh_patch_t* patch)
{
  // Y faces aren't communicated across z boundaries.
}

static void update_x1_zfaces(void* context, unimesh_t* mesh,
                             int i, int j, int k, real_t t,
                             field_metadata_t* md,
                             unimesh_patch_t* patch)
{
  // Z faces aren't communicated across x boundaries.
}

static void update_x2_zfaces(void* context, unimesh_t* mesh,
                             int i, int j, int k, real_t t,
                             field_metadata_t* md,
                             unimesh_patch_t* patch)
{
  // Z faces aren't communicated across x boundaries.
}

static void update_y1_zfaces(void* context, unimesh_t* mesh,
                             int i, int j, int k, real_t t,
                             field_metadata_t* md,
                             unimesh_patch_t* patch)
{
  // Z faces aren't communicated across y boundaries.
}

static void update_y2_zfaces(void* context, unimesh_t* mesh,
                             int i, int j, int k, real_t t,
                             field_metadata_t* md,
                             unimesh_patch_t* patch)
{
  // Z faces aren't communicated across y boundaries.
}

static void update_z1_zfaces(void* context, unimesh_t* mesh,
                             int i, int j, int k, real_t t,
                             field_metadata_t* md,
                             unimesh_patch_t* patch)
{
  constant_bc_t* bc = context;
  ASSERT(bc->num_components == patch->nc);

  DECLARE_UNIMESH_ZFACE_ARRAY(a, patch);
  for (int ii = 0; ii < patch->nx; ++ii)
    for (int jj = 0; jj < patch->ny; ++jj)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][jj][0][c] = bc->values[c];
}

static void update_z2_zfaces(void* context, unimesh_t* mesh,
                             int i, int j, int k, real_t t,
                             field_metadata_t* md,
                             unimesh_patch_t* patch)
{
  constant_bc_t* bc = context;
  ASSERT(bc->num_components == patch->nc);

  DECLARE_UNIMESH_ZFACE_ARRAY(a, patch);
  for (int ii = 0; ii < patch->nx; ++ii)
    for (int jj = 0; jj < patch->ny; ++jj)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][jj][patch->nz][c] = bc->values[c];
}

static void update_x1_xedges(void* context, unimesh_t* mesh,
                             int i, int j, int k, real_t t,
                             field_metadata_t* md,
                             unimesh_patch_t* patch)
{
  // X edges aren't communicated across x boundaries.
}

static void update_x2_xedges(void* context, unimesh_t* mesh,
                             int i, int j, int k, real_t t,
                             field_metadata_t* md,
                             unimesh_patch_t* patch)
{
  // X edges aren't communicated across x boundaries.
}

static void update_y1_xedges(void* context, unimesh_t* mesh,
                             int i, int j, int k, real_t t,
                             field_metadata_t* md,
                             unimesh_patch_t* patch)
{
  constant_bc_t* bc = context;
  ASSERT(bc->num_components == patch->nc);

  DECLARE_UNIMESH_XEDGE_ARRAY(a, patch);
  for (int ii = 0; ii < patch->nx; ++ii)
    for (int kk = 0; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][0][kk][c] = bc->values[c];
}

static void update_y2_xedges(void* context, unimesh_t* mesh,
                             int i, int j, int k, real_t t,
                             field_metadata_t* md,
                             unimesh_patch_t* patch)
{
  constant_bc_t* bc = context;
  ASSERT(bc->num_components == patch->nc);

  DECLARE_UNIMESH_XEDGE_ARRAY(a, patch);
  for (int ii = 0; ii < patch->nx; ++ii)
    for (int kk = 0; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][patch->ny][kk][c] = bc->values[c];
}

static void update_z1_xedges(void* context, unimesh_t* mesh,
                             int i, int j, int k, real_t t,
                             field_metadata_t* md,
                             unimesh_patch_t* patch)
{
  constant_bc_t* bc = context;
  ASSERT(bc->num_components == patch->nc);

  DECLARE_UNIMESH_XEDGE_ARRAY(a, patch);
  for (int ii = 0; ii < patch->nx; ++ii)
    for (int jj = 0; jj <= patch->ny; ++jj)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][jj][0][c] = bc->values[c];
}

static void update_z2_xedges(void* context, unimesh_t* mesh,
                             int i, int j, int k, real_t t,
                             field_metadata_t* md,
                             unimesh_patch_t* patch)
{
  constant_bc_t* bc = context;
  ASSERT(bc->num_components == patch->nc);

  DECLARE_UNIMESH_XEDGE_ARRAY(a, patch);
  for (int ii = 0; ii < patch->nx; ++ii)
    for (int jj = 0; jj <= patch->ny; ++jj)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][jj][patch->nz][c] = bc->values[c];
}

static void update_x1_yedges(void* context, unimesh_t* mesh,
                             int i, int j, int k, real_t t,
                             field_metadata_t* md,
                             unimesh_patch_t* patch)
{
  constant_bc_t* bc = context;
  ASSERT(bc->num_components == patch->nc);

  DECLARE_UNIMESH_YEDGE_ARRAY(a, patch);
  for (int jj = 0; jj < patch->ny; ++jj)
    for (int kk = 0; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[0][jj][kk][c] = bc->values[c];
}

static void update_x2_yedges(void* context, unimesh_t* mesh,
                             int i, int j, int k, real_t t,
                             field_metadata_t* md,
                             unimesh_patch_t* patch)
{
  constant_bc_t* bc = context;
  ASSERT(bc->num_components == patch->nc);

  DECLARE_UNIMESH_YEDGE_ARRAY(a, patch);
  for (int jj = 0; jj < patch->ny; ++jj)
    for (int kk = 0; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[patch->nx][jj][kk][c] = bc->values[c];
}

static void update_y1_yedges(void* context, unimesh_t* mesh,
                             int i, int j, int k, real_t t,
                             field_metadata_t* md,
                             unimesh_patch_t* patch)
{
  // Y edges aren't communicated across y boundaries.
}

static void update_y2_yedges(void* context, unimesh_t* mesh,
                             int i, int j, int k, real_t t,
                             field_metadata_t* md,
                             unimesh_patch_t* patch)
{
  // Y edges aren't communicated across y boundaries.
}

static void update_z1_yedges(void* context, unimesh_t* mesh,
                             int i, int j, int k, real_t t,
                             field_metadata_t* md,
                             unimesh_patch_t* patch)
{
  constant_bc_t* bc = context;
  ASSERT(bc->num_components == patch->nc);

  DECLARE_UNIMESH_YEDGE_ARRAY(a, patch);
  for (int ii = 0; ii <= patch->nx; ++ii)
    for (int jj = 0; jj < patch->ny; ++jj)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][jj][0][c] = bc->values[c];
}

static void update_z2_yedges(void* context, unimesh_t* mesh,
                             int i, int j, int k, real_t t,
                             field_metadata_t* md,
                             unimesh_patch_t* patch)
{
  constant_bc_t* bc = context;
  ASSERT(bc->num_components == patch->nc);

  DECLARE_UNIMESH_YEDGE_ARRAY(a, patch);
  for (int ii = 0; ii <= patch->nx; ++ii)
    for (int jj = 0; jj < patch->ny; ++jj)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][jj][patch->nz][c] = bc->values[c];
}

static void update_x1_zedges(void* context, unimesh_t* mesh,
                             int i, int j, int k, real_t t,
                             field_metadata_t* md,
                             unimesh_patch_t* patch)
{
  constant_bc_t* bc = context;
  ASSERT(bc->num_components == patch->nc);

  DECLARE_UNIMESH_ZEDGE_ARRAY(a, patch);
  for (int jj = 0; jj <= patch->ny; ++jj)
    for (int kk = 0; kk < patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[0][jj][kk][c] = bc->values[c];
}

static void update_x2_zedges(void* context, unimesh_t* mesh,
                             int i, int j, int k, real_t t,
                             field_metadata_t* md,
                             unimesh_patch_t* patch)
{
  constant_bc_t* bc = context;
  ASSERT(bc->num_components == patch->nc);

  DECLARE_UNIMESH_ZEDGE_ARRAY(a, patch);
  for (int jj = 0; jj <= patch->ny; ++jj)
    for (int kk = 0; kk < patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[patch->nx][jj][kk][c] = bc->values[c];
}

static void update_y1_zedges(void* context, unimesh_t* mesh,
                             int i, int j, int k, real_t t,
                             field_metadata_t* md,
                             unimesh_patch_t* patch)
{
  constant_bc_t* bc = context;
  ASSERT(bc->num_components == patch->nc);

  DECLARE_UNIMESH_ZEDGE_ARRAY(a, patch);
  for (int ii = 0; ii <= patch->nx; ++ii)
    for (int kk = 0; kk < patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][0][kk][c] = bc->values[c];
}

static void update_y2_zedges(void* context, unimesh_t* mesh,
                             int i, int j, int k, real_t t,
                             field_metadata_t* md,
                             unimesh_patch_t* patch)
{
  constant_bc_t* bc = context;
  ASSERT(bc->num_components == patch->nc);

  DECLARE_UNIMESH_ZEDGE_ARRAY(a, patch);
  for (int ii = 0; ii <= patch->nx; ++ii)
    for (int kk = 0; kk < patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][patch->ny][kk][c] = bc->values[c];
}

static void update_z1_zedges(void* context, unimesh_t* mesh,
                             int i, int j, int k, real_t t,
                             field_metadata_t* md,
                             unimesh_patch_t* patch)
{
  // Z edges aren't communicated across z boundaries.
}

static void update_z2_zedges(void* context, unimesh_t* mesh,
                             int i, int j, int k, real_t t,
                             field_metadata_t* md,
                             unimesh_patch_t* patch)
{
  // Z edges aren't communicated across z boundaries.
}

static void update_x1_nodes(void* context, unimesh_t* mesh,
                            int i, int j, int k, real_t t,
                            field_metadata_t* md,
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
                            field_metadata_t* md,
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
                            field_metadata_t* md,
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
                            field_metadata_t* md,
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
                            field_metadata_t* md,
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
                            field_metadata_t* md,
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
  vtable.start_update[1][0] = update_x1_xfaces;
  vtable.start_update[1][1] = update_x2_xfaces;
  vtable.start_update[1][2] = update_y1_xfaces;
  vtable.start_update[1][3] = update_y2_xfaces;
  vtable.start_update[1][4] = update_z1_xfaces;
  vtable.start_update[1][5] = update_z2_xfaces;
  vtable.start_update[2][0] = update_x1_yfaces;
  vtable.start_update[2][1] = update_x2_yfaces;
  vtable.start_update[2][2] = update_y1_yfaces;
  vtable.start_update[2][3] = update_y2_yfaces;
  vtable.start_update[2][4] = update_z1_yfaces;
  vtable.start_update[2][5] = update_z2_yfaces;
  vtable.start_update[3][0] = update_x1_zfaces;
  vtable.start_update[3][1] = update_x2_zfaces;
  vtable.start_update[3][2] = update_y1_zfaces;
  vtable.start_update[3][3] = update_y2_zfaces;
  vtable.start_update[3][4] = update_z1_zfaces;
  vtable.start_update[3][5] = update_z2_zfaces;
  vtable.start_update[4][0] = update_x1_xedges;
  vtable.start_update[4][1] = update_x2_xedges;
  vtable.start_update[4][2] = update_y1_xedges;
  vtable.start_update[4][3] = update_y2_xedges;
  vtable.start_update[4][4] = update_z1_xedges;
  vtable.start_update[4][5] = update_z2_xedges;
  vtable.start_update[5][0] = update_x1_yedges;
  vtable.start_update[5][1] = update_x2_yedges;
  vtable.start_update[5][2] = update_y1_yedges;
  vtable.start_update[5][3] = update_y2_yedges;
  vtable.start_update[5][4] = update_z1_yedges;
  vtable.start_update[5][5] = update_z2_yedges;
  vtable.start_update[6][0] = update_x1_zedges;
  vtable.start_update[6][1] = update_x2_zedges;
  vtable.start_update[6][2] = update_y1_zedges;
  vtable.start_update[6][3] = update_y2_zedges;
  vtable.start_update[6][4] = update_z1_zedges;
  vtable.start_update[6][5] = update_z2_zedges;
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
