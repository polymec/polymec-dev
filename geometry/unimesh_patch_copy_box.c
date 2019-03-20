// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "geometry/unimesh.h"
#include "geometry/unimesh_patch.h"

static void copy_cell(unimesh_patch_t* patch,
                      unimesh_patch_box_t* src_box,
                      unimesh_patch_box_t* dest_box,
                      unimesh_patch_t* dest)
{
  ASSERT(src_box->i2 <= patch->nx+2);
  ASSERT(src_box->j2 <= patch->ny+2);
  ASSERT(src_box->k2 <= patch->nz+2);
  ASSERT(dest_box->i2 <= dest->nx+2);
  ASSERT(dest_box->j2 <= dest->ny+2);
  ASSERT(dest_box->k2 <= dest->nz+2);

  DECLARE_UNIMESH_CELL_ARRAY(s, patch);
  DECLARE_UNIMESH_CELL_ARRAY(d, dest);

  int ni = src_box->i2 - src_box->i1;
  int nj = src_box->j2 - src_box->j1;
  int nk = src_box->k2 - src_box->k1;
  int si1 = src_box->i1, di1 = dest_box->i1;
  int sj1 = src_box->j1, dj1 = dest_box->j1;
  int sk1 = src_box->k1, dk1 = dest_box->k1;
  int nc = patch->nc;
  for (int i = 0; i < ni; ++i)
    for (int j = 0; j < nj; ++j)
      for (int k = 0; k < nk; ++k)
        for (int c = 0; c < nc; ++c)
          d[di1+i][dj1+j][dk1+k][c] = s[si1+i][sj1+j][sk1+k][c];
}

static void copy_xface(unimesh_patch_t* patch,
                       unimesh_patch_box_t* src_box,
                       unimesh_patch_box_t* dest_box,
                       unimesh_patch_t* dest)
{
  ASSERT(src_box->i2 <= patch->nx+1);
  ASSERT(src_box->j2 <= patch->ny);
  ASSERT(src_box->k2 <= patch->nz);
  ASSERT(dest_box->i2 <= dest->nx+1);
  ASSERT(dest_box->j2 <= dest->ny);
  ASSERT(dest_box->k2 <= dest->nz);

  DECLARE_UNIMESH_XFACE_ARRAY(s, patch);
  DECLARE_UNIMESH_XFACE_ARRAY(d, dest);

  int ni = src_box->i2 - src_box->i1;
  int nj = src_box->j2 - src_box->j1;
  int nk = src_box->k2 - src_box->k1;
  int si1 = src_box->i1, di1 = dest_box->i1;
  int sj1 = src_box->j1, dj1 = dest_box->j1;
  int sk1 = src_box->k1, dk1 = dest_box->k1;
  int nc = patch->nc;
  for (int i = 0; i < ni; ++i)
    for (int j = 0; j < nj; ++j)
      for (int k = 0; k < nk; ++k)
        for (int c = 0; c < nc; ++c)
          d[di1+i][dj1+j][dk1+k][c] = s[si1+i][sj1+j][sk1+k][c];
}

static void copy_yface(unimesh_patch_t* patch,
                       unimesh_patch_box_t* src_box,
                       unimesh_patch_box_t* dest_box,
                       unimesh_patch_t* dest)
{
  ASSERT(src_box->i2 <= patch->nx);
  ASSERT(src_box->j2 <= patch->ny+1);
  ASSERT(src_box->k2 <= patch->nz);
  ASSERT(dest_box->i2 <= dest->nx);
  ASSERT(dest_box->j2 <= dest->ny+1);
  ASSERT(dest_box->k2 <= dest->nz);

  DECLARE_UNIMESH_YFACE_ARRAY(s, patch);
  DECLARE_UNIMESH_YFACE_ARRAY(d, dest);

  int ni = src_box->i2 - src_box->i1;
  int nj = src_box->j2 - src_box->j1;
  int nk = src_box->k2 - src_box->k1;
  int si1 = src_box->i1, di1 = dest_box->i1;
  int sj1 = src_box->j1, dj1 = dest_box->j1;
  int sk1 = src_box->k1, dk1 = dest_box->k1;
  int nc = patch->nc;
  for (int i = 0; i < ni; ++i)
    for (int j = 0; j < nj; ++j)
      for (int k = 0; k < nk; ++k)
        for (int c = 0; c < nc; ++c)
          d[di1+i][dj1+j][dk1+k][c] = s[si1+i][sj1+j][sk1+k][c];
}

static void copy_zface(unimesh_patch_t* patch,
                       unimesh_patch_box_t* src_box,
                       unimesh_patch_box_t* dest_box,
                       unimesh_patch_t* dest)
{
  ASSERT(src_box->i2 <= patch->nx);
  ASSERT(src_box->j2 <= patch->ny);
  ASSERT(src_box->k2 <= patch->nz+1);
  ASSERT(dest_box->i2 <= dest->nx);
  ASSERT(dest_box->j2 <= dest->ny);
  ASSERT(dest_box->k2 <= dest->nz+1);

  DECLARE_UNIMESH_ZFACE_ARRAY(s, patch);
  DECLARE_UNIMESH_ZFACE_ARRAY(d, dest);

  int ni = src_box->i2 - src_box->i1;
  int nj = src_box->j2 - src_box->j1;
  int nk = src_box->k2 - src_box->k1;
  int si1 = src_box->i1, di1 = dest_box->i1;
  int sj1 = src_box->j1, dj1 = dest_box->j1;
  int sk1 = src_box->k1, dk1 = dest_box->k1;
  int nc = patch->nc;
  for (int i = 0; i < ni; ++i)
    for (int j = 0; j < nj; ++j)
      for (int k = 0; k < nk; ++k)
        for (int c = 0; c < nc; ++c)
          d[di1+i][dj1+j][dk1+k][c] = s[si1+i][sj1+j][sk1+k][c];
}

static void copy_xedge(unimesh_patch_t* patch,
                       unimesh_patch_box_t* src_box,
                       unimesh_patch_box_t* dest_box,
                       unimesh_patch_t* dest)
{
  ASSERT(src_box->i2 <= patch->nx);
  ASSERT(src_box->j2 <= patch->ny+1);
  ASSERT(src_box->k2 <= patch->nz+1);
  ASSERT(dest_box->i2 <= dest->nx);
  ASSERT(dest_box->j2 <= dest->ny+1);
  ASSERT(dest_box->k2 <= dest->nz+1);

  DECLARE_UNIMESH_XEDGE_ARRAY(s, patch);
  DECLARE_UNIMESH_XEDGE_ARRAY(d, dest);

  int ni = src_box->i2 - src_box->i1;
  int nj = src_box->j2 - src_box->j1;
  int nk = src_box->k2 - src_box->k1;
  int si1 = src_box->i1, di1 = dest_box->i1;
  int sj1 = src_box->j1, dj1 = dest_box->j1;
  int sk1 = src_box->k1, dk1 = dest_box->k1;
  int nc = patch->nc;
  for (int i = 0; i < ni; ++i)
    for (int j = 0; j < nj; ++j)
      for (int k = 0; k < nk; ++k)
        for (int c = 0; c < nc; ++c)
          d[di1+i][dj1+j][dk1+k][c] = s[si1+i][sj1+j][sk1+k][c];
}

static void copy_yedge(unimesh_patch_t* patch,
                       unimesh_patch_box_t* src_box,
                       unimesh_patch_box_t* dest_box,
                       unimesh_patch_t* dest)
{
  ASSERT(src_box->i2 <= patch->nx+1);
  ASSERT(src_box->j2 <= patch->ny);
  ASSERT(src_box->k2 <= patch->nz+1);
  ASSERT(dest_box->i2 <= dest->nx+1);
  ASSERT(dest_box->j2 <= dest->ny);
  ASSERT(dest_box->k2 <= dest->nz+1);

  DECLARE_UNIMESH_YEDGE_ARRAY(s, patch);
  DECLARE_UNIMESH_YEDGE_ARRAY(d, dest);

  int ni = src_box->i2 - src_box->i1;
  int nj = src_box->j2 - src_box->j1;
  int nk = src_box->k2 - src_box->k1;
  int si1 = src_box->i1, di1 = dest_box->i1;
  int sj1 = src_box->j1, dj1 = dest_box->j1;
  int sk1 = src_box->k1, dk1 = dest_box->k1;
  int nc = patch->nc;
  for (int i = 0; i < ni; ++i)
    for (int j = 0; j < nj; ++j)
      for (int k = 0; k < nk; ++k)
        for (int c = 0; c < nc; ++c)
          d[di1+i][dj1+j][dk1+k][c] = s[si1+i][sj1+j][sk1+k][c];
}

static void copy_zedge(unimesh_patch_t* patch,
                       unimesh_patch_box_t* src_box,
                       unimesh_patch_box_t* dest_box,
                       unimesh_patch_t* dest)
{
  ASSERT(src_box->i2 <= patch->nx+1);
  ASSERT(src_box->j2 <= patch->ny+1);
  ASSERT(src_box->k2 <= patch->nz);
  ASSERT(dest_box->i2 <= dest->nx+1);
  ASSERT(dest_box->j2 <= dest->ny+1);
  ASSERT(dest_box->k2 <= dest->nz);

  DECLARE_UNIMESH_ZEDGE_ARRAY(s, patch);
  DECLARE_UNIMESH_ZEDGE_ARRAY(d, dest);

  int ni = src_box->i2 - src_box->i1;
  int nj = src_box->j2 - src_box->j1;
  int nk = src_box->k2 - src_box->k1;
  int si1 = src_box->i1, di1 = dest_box->i1;
  int sj1 = src_box->j1, dj1 = dest_box->j1;
  int sk1 = src_box->k1, dk1 = dest_box->k1;
  int nc = patch->nc;
  for (int i = 0; i < ni; ++i)
    for (int j = 0; j < nj; ++j)
      for (int k = 0; k < nk; ++k)
        for (int c = 0; c < nc; ++c)
          d[di1+i][dj1+j][dk1+k][c] = s[si1+i][sj1+j][sk1+k][c];
}

static void copy_node(unimesh_patch_t* patch,
                      unimesh_patch_box_t* src_box,
                      unimesh_patch_box_t* dest_box,
                      unimesh_patch_t* dest)
{
  ASSERT(src_box->i2 <= patch->nx+1);
  ASSERT(src_box->j2 <= patch->ny+1);
  ASSERT(src_box->k2 <= patch->nz+1);
  ASSERT(dest_box->i2 <= dest->nx+1);
  ASSERT(dest_box->j2 <= dest->ny+1);
  ASSERT(dest_box->k2 <= dest->nz+1);

  DECLARE_UNIMESH_NODE_ARRAY(s, patch);
  DECLARE_UNIMESH_NODE_ARRAY(d, dest);

  int ni = src_box->i2 - src_box->i1;
  int nj = src_box->j2 - src_box->j1;
  int nk = src_box->k2 - src_box->k1;
  int si1 = src_box->i1, di1 = dest_box->i1;
  int sj1 = src_box->j1, dj1 = dest_box->j1;
  int sk1 = src_box->k1, dk1 = dest_box->k1;
  int nc = patch->nc;
  for (int i = 0; i < ni; ++i)
    for (int j = 0; j < nj; ++j)
      for (int k = 0; k < nk; ++k)
        for (int c = 0; c < nc; ++c)
          d[di1+i][dj1+j][dk1+k][c] = s[si1+i][sj1+j][sk1+k][c];
}

typedef void (*copy_box_func)(unimesh_patch_t* patch,
                              unimesh_patch_box_t* src_box,
                              unimesh_patch_box_t* dest_box,
                              unimesh_patch_t* dest);

void unimesh_patch_copy_box(unimesh_patch_t* patch,
                            unimesh_patch_box_t* src_box,
                            unimesh_patch_box_t* dest_box,
                            unimesh_patch_t* dest)
{
  ASSERT((src_box->i2 - src_box->i1) == (dest_box->i2 - dest_box->i1));
  ASSERT((src_box->j2 - src_box->j1) == (dest_box->j2 - dest_box->j1));
  ASSERT((src_box->k2 - src_box->k1) == (dest_box->k2 - dest_box->k1));
  ASSERT(src_box->i1 >= 0);
  ASSERT(src_box->j1 >= 0);
  ASSERT(src_box->k1 >= 0);
  ASSERT(dest_box->i1 >= 0);
  ASSERT(dest_box->j1 >= 0);
  ASSERT(dest_box->k1 >= 0);

  static copy_box_func copy[8] =
  {
    copy_cell, copy_xface, copy_yface, copy_zface,
    copy_xedge, copy_yedge, copy_zedge, copy_node
  };
  int c = (int)patch->centering;
  copy[c](patch, src_box, dest_box, dest);
}

