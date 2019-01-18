// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmocka.h"
#include "geometry/unimesh_patch.h"

static void test_unimesh_cell_patch(void** state) 
{ 
  unimesh_patch_t* patch = unimesh_patch_new(UNIMESH_CELL, 10, 10, 10, 4); 
  assert_true(patch->centering == UNIMESH_CELL);
  assert_int_equal(10, patch->nx);
  assert_int_equal(10, patch->ny);
  assert_int_equal(10, patch->nz);
  assert_int_equal(4, patch->nc);

  DECLARE_UNIMESH_CELL_ARRAY(a, patch);
  for (int i = 1; i <= patch->nx; ++i)
  {
    for (int j = 1; j <= patch->ny; ++j)
    {
      for (int k = 1; k <= patch->nz; ++k)
      {
        for (int l = 0; l < 4; ++l)
        {
          assert_true(reals_equal(a[i][j][k][l], 0.0));
          a[i][j][k][l] = (real_t)(12*12*4*i + 12*4*j + 4*k + l);
        }
      }
    }
  }

  unimesh_patch_box_t box;
  unimesh_patch_get_box(patch, &box);
  assert_int_equal(box.i1, 1);
  assert_int_equal(box.i2, patch->nx+1);
  assert_int_equal(box.j1, 1);
  assert_int_equal(box.j2, patch->ny+1);
  assert_int_equal(box.k1, 1);
  assert_int_equal(box.k2, patch->nz+1);

  unimesh_patch_get_boundary_box(patch, UNIMESH_X1_BOUNDARY, &box);
  assert_int_equal(box.i1, 0);
  assert_int_equal(box.i2, 1);
  unimesh_patch_get_boundary_box(patch, UNIMESH_X2_BOUNDARY, &box);
  assert_int_equal(box.i1, patch->nx+1);
  assert_int_equal(box.i2, patch->nx+2);
  unimesh_patch_get_boundary_box(patch, UNIMESH_Y1_BOUNDARY, &box);
  assert_int_equal(box.j1, 0);
  assert_int_equal(box.j2, 1);
  unimesh_patch_get_boundary_box(patch, UNIMESH_Y2_BOUNDARY, &box);
  assert_int_equal(box.j1, patch->ny+1);
  assert_int_equal(box.j2, patch->ny+2);
  unimesh_patch_get_boundary_box(patch, UNIMESH_Z1_BOUNDARY, &box);
  assert_int_equal(box.k1, 0);
  assert_int_equal(box.k2, 1);
  unimesh_patch_get_boundary_box(patch, UNIMESH_Z2_BOUNDARY, &box);
  assert_int_equal(box.k1, patch->nz+1);
  assert_int_equal(box.k2, patch->nz+2);

  unimesh_patch_t* patch1 = unimesh_patch_clone(patch);
  assert_true(patch1->centering == UNIMESH_CELL);
  assert_int_equal(10, patch1->nx);
  assert_int_equal(10, patch1->ny);
  assert_int_equal(10, patch1->nz);
  assert_int_equal(4, patch1->nc);

  DECLARE_UNIMESH_CELL_ARRAY(a1, patch);
  for (int i = 1; i <= patch->nx; ++i)
    for (int j = 1; j <= patch->ny; ++j)
      for (int k = 1; k <= patch->nz; ++k)
        for (int l = 0; l < 4; ++l)
          assert_true(reals_equal(a1[i][j][k][l], a[i][j][k][l]));

  unimesh_patch_free(patch1); 
  unimesh_patch_free(patch); 
} 

static void test_unimesh_xface_patch(void** state) 
{ 
  unimesh_patch_t* patch = unimesh_patch_new(UNIMESH_XFACE, 10, 10, 10, 4); 
  assert_true(patch->centering == UNIMESH_XFACE);
  assert_int_equal(10, patch->nx);
  assert_int_equal(10, patch->ny);
  assert_int_equal(10, patch->nz);
  assert_int_equal(4, patch->nc);

  DECLARE_UNIMESH_XFACE_ARRAY(a, patch);
  for (int i = 0; i <= patch->nx; ++i)
  {
    for (int j = 0; j < patch->ny; ++j)
    {
      for (int k = 0; k < patch->nz; ++k)
      {
        for (int l = 0; l < 4; ++l)
        {
          assert_true(reals_equal(a[i][j][k][l], 0.0));
          a[i][j][k][l] = (real_t)(12*12*4*i + 12*4*j + 4*k + l);
        }
      }
    }
  }

  unimesh_patch_box_t box;
  unimesh_patch_get_box(patch, &box);
  assert_int_equal(box.i1, 0);
  assert_int_equal(box.i2, patch->nx+1);
  assert_int_equal(box.j1, 0);
  assert_int_equal(box.j2, patch->ny);
  assert_int_equal(box.k1, 0);
  assert_int_equal(box.k2, patch->nz);

  unimesh_patch_get_boundary_box(patch, UNIMESH_X1_BOUNDARY, &box);
  assert_int_equal(box.i1, 0);
  assert_int_equal(box.i2, 1);
  unimesh_patch_get_boundary_box(patch, UNIMESH_X2_BOUNDARY, &box);
  assert_int_equal(box.i1, patch->nx);
  assert_int_equal(box.i2, patch->nx+1);
  unimesh_patch_get_boundary_box(patch, UNIMESH_Y1_BOUNDARY, &box);
  assert_int_equal(box.j1, 0);
  assert_int_equal(box.j2, 1);
  unimesh_patch_get_boundary_box(patch, UNIMESH_Y2_BOUNDARY, &box);
  assert_int_equal(box.j1, patch->ny-1);
  assert_int_equal(box.j2, patch->ny);
  unimesh_patch_get_boundary_box(patch, UNIMESH_Z1_BOUNDARY, &box);
  assert_int_equal(box.k1, 0);
  assert_int_equal(box.k2, 1);
  unimesh_patch_get_boundary_box(patch, UNIMESH_Z2_BOUNDARY, &box);
  assert_int_equal(box.k1, patch->nz-1);
  assert_int_equal(box.k2, patch->nz);

  unimesh_patch_t* patch1 = unimesh_patch_clone(patch);
  assert_true(patch1->centering == UNIMESH_XFACE);
  assert_int_equal(10, patch1->nx);
  assert_int_equal(10, patch1->ny);
  assert_int_equal(10, patch1->nz);
  assert_int_equal(4, patch1->nc);

  DECLARE_UNIMESH_XFACE_ARRAY(a1, patch);
  for (int i = 0; i <= patch->nx; ++i)
    for (int j = 0; j < patch->ny; ++j)
      for (int k = 0; k < patch->nz; ++k)
        for (int l = 0; l < 4; ++l)
          assert_true(reals_equal(a1[i][j][k][l], a[i][j][k][l]));

  unimesh_patch_free(patch1); 
  unimesh_patch_free(patch); 
} 

static void test_unimesh_yface_patch(void** state) 
{ 
  unimesh_patch_t* patch = unimesh_patch_new(UNIMESH_YFACE, 10, 10, 10, 4); 
  assert_true(patch->centering == UNIMESH_YFACE);
  assert_int_equal(10, patch->nx);
  assert_int_equal(10, patch->ny);
  assert_int_equal(10, patch->nz);
  assert_int_equal(4, patch->nc);

  DECLARE_UNIMESH_YFACE_ARRAY(a, patch);
  for (int i = 0; i < patch->nx; ++i)
  {
    for (int j = 0; j <= patch->ny; ++j)
    {
      for (int k = 0; k < patch->nz; ++k)
      {
        for (int l = 0; l < 4; ++l)
        {
          assert_true(reals_equal(a[i][j][k][l], 0.0));
          a[i][j][k][l] = (real_t)(12*12*4*i + 12*4*j + 4*k + l);
        }
      }
    }
  }

  unimesh_patch_box_t box;
  unimesh_patch_get_box(patch, &box);
  assert_int_equal(box.i1, 0);
  assert_int_equal(box.i2, patch->nx);
  assert_int_equal(box.j1, 0);
  assert_int_equal(box.j2, patch->ny+1);
  assert_int_equal(box.k1, 0);
  assert_int_equal(box.k2, patch->nz);

  unimesh_patch_get_boundary_box(patch, UNIMESH_X1_BOUNDARY, &box);
  assert_int_equal(box.i1, 0);
  assert_int_equal(box.i2, 1);
  unimesh_patch_get_boundary_box(patch, UNIMESH_X2_BOUNDARY, &box);
  assert_int_equal(box.i1, patch->nx-1);
  assert_int_equal(box.i2, patch->nx);
  unimesh_patch_get_boundary_box(patch, UNIMESH_Y1_BOUNDARY, &box);
  assert_int_equal(box.j1, 0);
  assert_int_equal(box.j2, 1);
  unimesh_patch_get_boundary_box(patch, UNIMESH_Y2_BOUNDARY, &box);
  assert_int_equal(box.j1, patch->ny);
  assert_int_equal(box.j2, patch->ny+1);
  unimesh_patch_get_boundary_box(patch, UNIMESH_Z1_BOUNDARY, &box);
  assert_int_equal(box.k1, 0);
  assert_int_equal(box.k2, 1);
  unimesh_patch_get_boundary_box(patch, UNIMESH_Z2_BOUNDARY, &box);
  assert_int_equal(box.k1, patch->nz-1);
  assert_int_equal(box.k2, patch->nz);

  unimesh_patch_t* patch1 = unimesh_patch_clone(patch);
  assert_true(patch1->centering == UNIMESH_YFACE);
  assert_int_equal(10, patch1->nx);
  assert_int_equal(10, patch1->ny);
  assert_int_equal(10, patch1->nz);
  assert_int_equal(4, patch1->nc);

  DECLARE_UNIMESH_YFACE_ARRAY(a1, patch);
  for (int i = 0; i < patch->nx; ++i)
    for (int j = 0; j <= patch->ny; ++j)
      for (int k = 0; k < patch->nz; ++k)
        for (int l = 0; l < 4; ++l)
          assert_true(reals_equal(a1[i][j][k][l], a[i][j][k][l]));

  unimesh_patch_free(patch1); 
  unimesh_patch_free(patch); 
} 

static void test_unimesh_zface_patch(void** state) 
{ 
  unimesh_patch_t* patch = unimesh_patch_new(UNIMESH_ZFACE, 10, 10, 10, 4); 
  assert_true(patch->centering == UNIMESH_ZFACE);
  assert_int_equal(10, patch->nx);
  assert_int_equal(10, patch->ny);
  assert_int_equal(10, patch->nz);
  assert_int_equal(4, patch->nc);

  DECLARE_UNIMESH_ZFACE_ARRAY(a, patch);
  for (int i = 0; i < patch->nx; ++i)
  {
    for (int j = 0; j < patch->ny; ++j)
    {
      for (int k = 0; k <= patch->nz; ++k)
      {
        for (int l = 0; l < 4; ++l)
        {
          assert_true(reals_equal(a[i][j][k][l], 0.0));
          a[i][j][k][l] = (real_t)(12*12*4*i + 12*4*j + 4*k + l);
        }
      }
    }
  }

  unimesh_patch_box_t box;
  unimesh_patch_get_box(patch, &box);
  assert_int_equal(box.i1, 0);
  assert_int_equal(box.i2, patch->nx);
  assert_int_equal(box.j1, 0);
  assert_int_equal(box.j2, patch->ny);
  assert_int_equal(box.k1, 0);
  assert_int_equal(box.k2, patch->nz+1);

  unimesh_patch_get_boundary_box(patch, UNIMESH_X1_BOUNDARY, &box);
  assert_int_equal(box.i1, 0);
  assert_int_equal(box.i2, 1);
  unimesh_patch_get_boundary_box(patch, UNIMESH_X2_BOUNDARY, &box);
  assert_int_equal(box.i1, patch->nx-1);
  assert_int_equal(box.i2, patch->nx);
  unimesh_patch_get_boundary_box(patch, UNIMESH_Y1_BOUNDARY, &box);
  assert_int_equal(box.j1, 0);
  assert_int_equal(box.j2, 1);
  unimesh_patch_get_boundary_box(patch, UNIMESH_Y2_BOUNDARY, &box);
  assert_int_equal(box.j1, patch->ny-1);
  assert_int_equal(box.j2, patch->ny);
  unimesh_patch_get_boundary_box(patch, UNIMESH_Z1_BOUNDARY, &box);
  assert_int_equal(box.k1, 0);
  assert_int_equal(box.k2, 1);
  unimesh_patch_get_boundary_box(patch, UNIMESH_Z2_BOUNDARY, &box);
  assert_int_equal(box.k1, patch->nz);
  assert_int_equal(box.k2, patch->nz+1);

  unimesh_patch_t* patch1 = unimesh_patch_clone(patch);
  assert_true(patch1->centering == UNIMESH_ZFACE);
  assert_int_equal(10, patch1->nx);
  assert_int_equal(10, patch1->ny);
  assert_int_equal(10, patch1->nz);
  assert_int_equal(4, patch1->nc);

  DECLARE_UNIMESH_ZFACE_ARRAY(a1, patch);
  for (int i = 0; i < patch->nx; ++i)
    for (int j = 0; j < patch->ny; ++j)
      for (int k = 0; k <= patch->nz; ++k)
        for (int l = 0; l < 4; ++l)
          assert_true(reals_equal(a1[i][j][k][l], a[i][j][k][l]));

  unimesh_patch_free(patch1); 
  unimesh_patch_free(patch); 
} 

static void test_unimesh_xedge_patch(void** state) 
{ 
  unimesh_patch_t* patch = unimesh_patch_new(UNIMESH_XEDGE, 10, 10, 10, 4); 
  assert_true(patch->centering == UNIMESH_XEDGE);
  assert_int_equal(10, patch->nx);
  assert_int_equal(10, patch->ny);
  assert_int_equal(10, patch->nz);
  assert_int_equal(4, patch->nc);

  DECLARE_UNIMESH_XEDGE_ARRAY(a, patch);
  for (int i = 0; i < patch->nx; ++i)
  {
    for (int j = 0; j <= patch->ny; ++j)
    {
      for (int k = 0; k <= patch->nz; ++k)
      {
        for (int l = 0; l < 4; ++l)
        {
          assert_true(reals_equal(a[i][j][k][l], 0.0));
          a[i][j][k][l] = (real_t)(12*12*4*i + 12*4*j + 4*k + l);
        }
      }
    }
  }

  unimesh_patch_box_t box;
  unimesh_patch_get_box(patch, &box);
  assert_int_equal(box.i1, 0);
  assert_int_equal(box.i2, patch->nx);
  assert_int_equal(box.j1, 0);
  assert_int_equal(box.j2, patch->ny+1);
  assert_int_equal(box.k1, 0);
  assert_int_equal(box.k2, patch->nz+1);

  unimesh_patch_get_boundary_box(patch, UNIMESH_X1_BOUNDARY, &box);
  assert_int_equal(box.i1, 0);
  assert_int_equal(box.i2, 1);
  unimesh_patch_get_boundary_box(patch, UNIMESH_X2_BOUNDARY, &box);
  assert_int_equal(box.i1, patch->nx-1);
  assert_int_equal(box.i2, patch->nx);
  unimesh_patch_get_boundary_box(patch, UNIMESH_Y1_BOUNDARY, &box);
  assert_int_equal(box.j1, 0);
  assert_int_equal(box.j2, 1);
  unimesh_patch_get_boundary_box(patch, UNIMESH_Y2_BOUNDARY, &box);
  assert_int_equal(box.j1, patch->ny);
  assert_int_equal(box.j2, patch->ny+1);
  unimesh_patch_get_boundary_box(patch, UNIMESH_Z1_BOUNDARY, &box);
  assert_int_equal(box.k1, 0);
  assert_int_equal(box.k2, 1);
  unimesh_patch_get_boundary_box(patch, UNIMESH_Z2_BOUNDARY, &box);
  assert_int_equal(box.k1, patch->nz);
  assert_int_equal(box.k2, patch->nz+1);

  unimesh_patch_t* patch1 = unimesh_patch_clone(patch);
  assert_true(patch1->centering == UNIMESH_XEDGE);
  assert_int_equal(10, patch1->nx);
  assert_int_equal(10, patch1->ny);
  assert_int_equal(10, patch1->nz);
  assert_int_equal(4, patch1->nc);

  DECLARE_UNIMESH_XEDGE_ARRAY(a1, patch);
  for (int i = 0; i < patch->nx; ++i)
    for (int j = 0; j <= patch->ny; ++j)
      for (int k = 0; k <= patch->nz; ++k)
        for (int l = 0; l < 4; ++l)
          assert_true(reals_equal(a1[i][j][k][l], a[i][j][k][l]));

  unimesh_patch_free(patch1); 
  unimesh_patch_free(patch); 
} 

static void test_unimesh_yedge_patch(void** state) 
{ 
  unimesh_patch_t* patch = unimesh_patch_new(UNIMESH_YEDGE, 10, 10, 10, 4); 
  assert_true(patch->centering == UNIMESH_YEDGE);
  assert_int_equal(10, patch->nx);
  assert_int_equal(10, patch->ny);
  assert_int_equal(10, patch->nz);
  assert_int_equal(4, patch->nc);

  DECLARE_UNIMESH_YEDGE_ARRAY(a, patch);
  for (int i = 0; i <= patch->nx; ++i)
  {
    for (int j = 0; j < patch->ny; ++j)
    {
      for (int k = 0; k <= patch->nz; ++k)
      {
        for (int l = 0; l < 4; ++l)
        {
          assert_true(reals_equal(a[i][j][k][l], 0.0));
          a[i][j][k][l] = (real_t)(12*12*4*i + 12*4*j + 4*k + l);
        }
      }
    }
  }

  unimesh_patch_box_t box;
  unimesh_patch_get_box(patch, &box);
  assert_int_equal(box.i1, 0);
  assert_int_equal(box.i2, patch->nx+1);
  assert_int_equal(box.j1, 0);
  assert_int_equal(box.j2, patch->ny);
  assert_int_equal(box.k1, 0);
  assert_int_equal(box.k2, patch->nz+1);

  unimesh_patch_get_boundary_box(patch, UNIMESH_X1_BOUNDARY, &box);
  assert_int_equal(box.i1, 0);
  assert_int_equal(box.i2, 1);
  unimesh_patch_get_boundary_box(patch, UNIMESH_X2_BOUNDARY, &box);
  assert_int_equal(box.i1, patch->nx);
  assert_int_equal(box.i2, patch->nx+1);
  unimesh_patch_get_boundary_box(patch, UNIMESH_Y1_BOUNDARY, &box);
  assert_int_equal(box.j1, 0);
  assert_int_equal(box.j2, 1);
  unimesh_patch_get_boundary_box(patch, UNIMESH_Y2_BOUNDARY, &box);
  assert_int_equal(box.j1, patch->ny-1);
  assert_int_equal(box.j2, patch->ny);
  unimesh_patch_get_boundary_box(patch, UNIMESH_Z1_BOUNDARY, &box);
  assert_int_equal(box.k1, 0);
  assert_int_equal(box.k2, 1);
  unimesh_patch_get_boundary_box(patch, UNIMESH_Z2_BOUNDARY, &box);
  assert_int_equal(box.k1, patch->nz);
  assert_int_equal(box.k2, patch->nz+1);

  unimesh_patch_t* patch1 = unimesh_patch_clone(patch);
  assert_true(patch1->centering == UNIMESH_YEDGE);
  assert_int_equal(10, patch1->nx);
  assert_int_equal(10, patch1->ny);
  assert_int_equal(10, patch1->nz);
  assert_int_equal(4, patch1->nc);

  DECLARE_UNIMESH_YEDGE_ARRAY(a1, patch);
  for (int i = 0; i <= patch->nx; ++i)
    for (int j = 0; j < patch->ny; ++j)
      for (int k = 0; k <= patch->nz; ++k)
        for (int l = 0; l < 4; ++l)
          assert_true(reals_equal(a1[i][j][k][l], a[i][j][k][l]));

  unimesh_patch_free(patch1); 
  unimesh_patch_free(patch); 
} 

static void test_unimesh_zedge_patch(void** state) 
{ 
  unimesh_patch_t* patch = unimesh_patch_new(UNIMESH_ZEDGE, 10, 10, 10, 4); 
  assert_true(patch->centering == UNIMESH_ZEDGE);
  assert_int_equal(10, patch->nx);
  assert_int_equal(10, patch->ny);
  assert_int_equal(10, patch->nz);
  assert_int_equal(4, patch->nc);

  DECLARE_UNIMESH_ZEDGE_ARRAY(a, patch);
  for (int i = 0; i <= patch->nx; ++i)
  {
    for (int j = 0; j <= patch->ny; ++j)
    {
      for (int k = 0; k < patch->nz; ++k)
      {
        for (int l = 0; l < 4; ++l)
        {
          assert_true(reals_equal(a[i][j][k][l], 0.0));
          a[i][j][k][l] = (real_t)(12*12*4*i + 12*4*j + 4*k + l);
        }
      }
    }
  }

  unimesh_patch_box_t box;
  unimesh_patch_get_box(patch, &box);
  assert_int_equal(box.i1, 0);
  assert_int_equal(box.i2, patch->nx+1);
  assert_int_equal(box.j1, 0);
  assert_int_equal(box.j2, patch->ny+1);
  assert_int_equal(box.k1, 0);
  assert_int_equal(box.k2, patch->nz);

  unimesh_patch_get_boundary_box(patch, UNIMESH_X1_BOUNDARY, &box);
  assert_int_equal(box.i1, 0);
  assert_int_equal(box.i2, 1);
  unimesh_patch_get_boundary_box(patch, UNIMESH_X2_BOUNDARY, &box);
  assert_int_equal(box.i1, patch->nx);
  assert_int_equal(box.i2, patch->nx+1);
  unimesh_patch_get_boundary_box(patch, UNIMESH_Y1_BOUNDARY, &box);
  assert_int_equal(box.j1, 0);
  assert_int_equal(box.j2, 1);
  unimesh_patch_get_boundary_box(patch, UNIMESH_Y2_BOUNDARY, &box);
  assert_int_equal(box.j1, patch->ny);
  assert_int_equal(box.j2, patch->ny+1);
  unimesh_patch_get_boundary_box(patch, UNIMESH_Z1_BOUNDARY, &box);
  assert_int_equal(box.k1, 0);
  assert_int_equal(box.k2, 1);
  unimesh_patch_get_boundary_box(patch, UNIMESH_Z2_BOUNDARY, &box);
  assert_int_equal(box.k1, patch->nz-1);
  assert_int_equal(box.k2, patch->nz);

  unimesh_patch_t* patch1 = unimesh_patch_clone(patch);
  assert_true(patch1->centering == UNIMESH_ZEDGE);
  assert_int_equal(10, patch1->nx);
  assert_int_equal(10, patch1->ny);
  assert_int_equal(10, patch1->nz);
  assert_int_equal(4, patch1->nc);

  DECLARE_UNIMESH_ZEDGE_ARRAY(a1, patch);
  for (int i = 0; i <= patch->nx; ++i)
    for (int j = 0; j <= patch->ny; ++j)
      for (int k = 0; k < patch->nz; ++k)
        for (int l = 0; l < 4; ++l)
          assert_true(reals_equal(a1[i][j][k][l], a[i][j][k][l]));

  unimesh_patch_free(patch1); 
  unimesh_patch_free(patch); 
} 

static void test_unimesh_node_patch(void** state) 
{ 
  unimesh_patch_t* patch = unimesh_patch_new(UNIMESH_NODE, 10, 10, 10, 4); 
  assert_true(patch->centering == UNIMESH_NODE);
  assert_int_equal(10, patch->nx);
  assert_int_equal(10, patch->ny);
  assert_int_equal(10, patch->nz);
  assert_int_equal(4, patch->nc);

  DECLARE_UNIMESH_NODE_ARRAY(a, patch);
  for (int i = 0; i <= patch->nx; ++i)
  {
    for (int j = 0; j <= patch->ny; ++j)
    {
      for (int k = 0; k <= patch->nz; ++k)
      {
        for (int l = 0; l < 4; ++l)
        {
          assert_true(reals_equal(a[i][j][k][l], 0.0));
          a[i][j][k][l] = (real_t)(12*12*4*i + 12*4*j + 4*k + l);
        }
      }
    }
  }

  unimesh_patch_box_t box;
  unimesh_patch_get_box(patch, &box);
  assert_int_equal(box.i1, 0);
  assert_int_equal(box.i2, patch->nx+1);
  assert_int_equal(box.j1, 0);
  assert_int_equal(box.j2, patch->ny+1);
  assert_int_equal(box.k1, 0);
  assert_int_equal(box.k2, patch->nz+1);

  unimesh_patch_get_boundary_box(patch, UNIMESH_X1_BOUNDARY, &box);
  assert_int_equal(box.i1, 0);
  assert_int_equal(box.i2, 1);
  unimesh_patch_get_boundary_box(patch, UNIMESH_X2_BOUNDARY, &box);
  assert_int_equal(box.i1, patch->nx);
  assert_int_equal(box.i2, patch->nx+1);
  unimesh_patch_get_boundary_box(patch, UNIMESH_Y1_BOUNDARY, &box);
  assert_int_equal(box.j1, 0);
  assert_int_equal(box.j2, 1);
  unimesh_patch_get_boundary_box(patch, UNIMESH_Y2_BOUNDARY, &box);
  assert_int_equal(box.j1, patch->ny);
  assert_int_equal(box.j2, patch->ny+1);
  unimesh_patch_get_boundary_box(patch, UNIMESH_Z1_BOUNDARY, &box);
  assert_int_equal(box.k1, 0);
  assert_int_equal(box.k2, 1);
  unimesh_patch_get_boundary_box(patch, UNIMESH_Z2_BOUNDARY, &box);
  assert_int_equal(box.k1, patch->nz);
  assert_int_equal(box.k2, patch->nz+1);

  unimesh_patch_t* patch1 = unimesh_patch_clone(patch);
  assert_true(patch1->centering == UNIMESH_NODE);
  assert_int_equal(10, patch1->nx);
  assert_int_equal(10, patch1->ny);
  assert_int_equal(10, patch1->nz);
  assert_int_equal(4, patch1->nc);

  DECLARE_UNIMESH_NODE_ARRAY(a1, patch);
  for (int i = 0; i <= patch->nx; ++i)
    for (int j = 0; j <= patch->ny; ++j)
      for (int k = 0; k <= patch->nz; ++k)
        for (int l = 0; l < 4; ++l)
          assert_true(reals_equal(a1[i][j][k][l], a[i][j][k][l]));

  unimesh_patch_free(patch1); 
  unimesh_patch_free(patch); 
} 

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_unimesh_cell_patch),
    cmocka_unit_test(test_unimesh_xface_patch),
    cmocka_unit_test(test_unimesh_yface_patch),
    cmocka_unit_test(test_unimesh_zface_patch),
    cmocka_unit_test(test_unimesh_xedge_patch),
    cmocka_unit_test(test_unimesh_yedge_patch),
    cmocka_unit_test(test_unimesh_zedge_patch),
    cmocka_unit_test(test_unimesh_node_patch)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
