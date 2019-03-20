// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "geometry/unimesh.h"
#include "geometry/unimesh_patch.h"

static void fill_x1_cell(unimesh_patch_t* patch, real_t* data)
{
  DECLARE_UNIMESH_CELL_ARRAY(a, patch);
  for (int jj = 1; jj <= patch->ny; ++jj)
    for (int kk = 1; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[0][jj][kk][c] = data[c];
}

static void fill_x2_cell(unimesh_patch_t* patch, real_t* data)
{
  DECLARE_UNIMESH_CELL_ARRAY(a, patch);
  for (int jj = 1; jj <= patch->ny; ++jj)
    for (int kk = 1; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[patch->nx+1][jj][kk][c] = data[c];
}

static void fill_y1_cell(unimesh_patch_t* patch, real_t* data)
{
  DECLARE_UNIMESH_CELL_ARRAY(a, patch);
  for (int ii = 1; ii <= patch->nx; ++ii)
    for (int kk = 1; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][0][kk][c] = data[c];
}

static void fill_y2_cell(unimesh_patch_t* patch, real_t* data)
{
  DECLARE_UNIMESH_CELL_ARRAY(a, patch);
  for (int ii = 1; ii <= patch->nx; ++ii)
    for (int kk = 1; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][patch->ny+1][kk][c] = data[c];
}

static void fill_z1_cell(unimesh_patch_t* patch, real_t* data)
{
  DECLARE_UNIMESH_CELL_ARRAY(a, patch);
  for (int ii = 1; ii <= patch->nx; ++ii)
    for (int jj = 1; jj <= patch->ny; ++jj)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][jj][0][c] = data[c];
}

static void fill_z2_cell(unimesh_patch_t* patch, real_t* data)
{
  DECLARE_UNIMESH_CELL_ARRAY(a, patch);
  for (int ii = 1; ii <= patch->nx; ++ii)
    for (int jj = 1; jj <= patch->ny; ++jj)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][jj][patch->nz+1][c] = data[c];
}

static void fill_x1_xface(unimesh_patch_t* patch, real_t* data)
{
  DECLARE_UNIMESH_XFACE_ARRAY(a, patch);
  for (int jj = 0; jj < patch->ny; ++jj)
    for (int kk = 0; kk < patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[0][jj][kk][c] = data[c];
}

static void fill_x2_xface(unimesh_patch_t* patch, real_t* data)
{
  DECLARE_UNIMESH_XFACE_ARRAY(a, patch);
  for (int jj = 0; jj < patch->ny; ++jj)
    for (int kk = 0; kk < patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[patch->nx][jj][kk][c] = data[c];
}

static void fill_y1_xface(unimesh_patch_t* patch, real_t* data)
{
}

static void fill_y2_xface(unimesh_patch_t* patch, real_t* data)
{
}

static void fill_z1_xface(unimesh_patch_t* patch, real_t* data)
{
}

static void fill_z2_xface(unimesh_patch_t* patch, real_t* data)
{
}

static void fill_x1_yface(unimesh_patch_t* patch, real_t* data)
{
}

static void fill_x2_yface(unimesh_patch_t* patch, real_t* data)
{
}

static void fill_y1_yface(unimesh_patch_t* patch, real_t* data)
{
  DECLARE_UNIMESH_YFACE_ARRAY(a, patch);
  for (int ii = 0; ii < patch->nx; ++ii)
    for (int kk = 0; kk < patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][0][kk][c] = data[c];
}

static void fill_y2_yface(unimesh_patch_t* patch, real_t* data)
{
  DECLARE_UNIMESH_YFACE_ARRAY(a, patch);
  for (int ii = 0; ii < patch->nx; ++ii)
    for (int kk = 0; kk < patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][patch->ny][kk][c] = data[c];
}

static void fill_z1_yface(unimesh_patch_t* patch, real_t* data)
{
}

static void fill_z2_yface(unimesh_patch_t* patch, real_t* data)
{
}

static void fill_x1_zface(unimesh_patch_t* patch, real_t* data)
{
}

static void fill_x2_zface(unimesh_patch_t* patch, real_t* data)
{
}

static void fill_y1_zface(unimesh_patch_t* patch, real_t* data)
{
}

static void fill_y2_zface(unimesh_patch_t* patch, real_t* data)
{
}

static void fill_z1_zface(unimesh_patch_t* patch, real_t* data)
{
  DECLARE_UNIMESH_ZFACE_ARRAY(a, patch);
  for (int ii = 0; ii < patch->nx; ++ii)
    for (int jj = 0; jj < patch->ny; ++jj)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][jj][0][c] = data[c];
}

static void fill_z2_zface(unimesh_patch_t* patch, real_t* data)
{
  DECLARE_UNIMESH_ZFACE_ARRAY(a, patch);
  for (int ii = 0; ii < patch->nx; ++ii)
    for (int jj = 0; jj < patch->ny; ++jj)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][jj][patch->nz][c] = data[c];
}

static void fill_x1_xedge(unimesh_patch_t* patch, real_t* data)
{
}

static void fill_x2_xedge(unimesh_patch_t* patch, real_t* data)
{
}

static void fill_y1_xedge(unimesh_patch_t* patch, real_t* data)
{
  DECLARE_UNIMESH_XEDGE_ARRAY(a, patch);
  for (int ii = 0; ii < patch->nx; ++ii)
    for (int kk = 0; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][0][kk][c] = data[c];
}

static void fill_y2_xedge(unimesh_patch_t* patch, real_t* data)
{
  DECLARE_UNIMESH_XEDGE_ARRAY(a, patch);
  for (int ii = 0; ii < patch->nx; ++ii)
    for (int kk = 0; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][patch->ny][kk][c] = data[c];
}

static void fill_z1_xedge(unimesh_patch_t* patch, real_t* data)
{
  DECLARE_UNIMESH_XEDGE_ARRAY(a, patch);
  for (int ii = 0; ii < patch->nx; ++ii)
    for (int jj = 0; jj <= patch->ny; ++jj)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][jj][0][c] = data[c];
}

static void fill_z2_xedge(unimesh_patch_t* patch, real_t* data)
{
  DECLARE_UNIMESH_XEDGE_ARRAY(a, patch);
  for (int ii = 0; ii < patch->nx; ++ii)
    for (int jj = 0; jj <= patch->ny; ++jj)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][jj][patch->nz][c] = data[c];
}

static void fill_x1_yedge(unimesh_patch_t* patch, real_t* data)
{
  DECLARE_UNIMESH_YEDGE_ARRAY(a, patch);
  for (int jj = 0; jj < patch->ny; ++jj)
    for (int kk = 0; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[0][jj][kk][c] = data[c];
}

static void fill_x2_yedge(unimesh_patch_t* patch, real_t* data)
{
  DECLARE_UNIMESH_YEDGE_ARRAY(a, patch);
  for (int jj = 0; jj < patch->ny; ++jj)
    for (int kk = 0; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[patch->nx][jj][kk][c] = data[c];
}

static void fill_y1_yedge(unimesh_patch_t* patch, real_t* data)
{
}

static void fill_y2_yedge(unimesh_patch_t* patch, real_t* data)
{
}

static void fill_z1_yedge(unimesh_patch_t* patch, real_t* data)
{
  DECLARE_UNIMESH_YEDGE_ARRAY(a, patch);
  for (int ii = 0; ii <= patch->nx; ++ii)
    for (int jj = 0; jj < patch->ny; ++jj)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][jj][0][c] = data[c];
}

static void fill_z2_yedge(unimesh_patch_t* patch, real_t* data)
{
  DECLARE_UNIMESH_YEDGE_ARRAY(a, patch);
  for (int ii = 0; ii <= patch->nx; ++ii)
    for (int jj = 0; jj < patch->ny; ++jj)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][jj][patch->nz][c] = data[c];
}

static void fill_x1_zedge(unimesh_patch_t* patch, real_t* data)
{
  DECLARE_UNIMESH_ZEDGE_ARRAY(a, patch);
  for (int jj = 0; jj <= patch->ny; ++jj)
    for (int kk = 0; kk < patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[0][jj][kk][c] = data[c];
}

static void fill_x2_zedge(unimesh_patch_t* patch, real_t* data)
{
  DECLARE_UNIMESH_ZEDGE_ARRAY(a, patch);
  for (int jj = 0; jj <= patch->ny; ++jj)
    for (int kk = 0; kk < patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[patch->nx][jj][kk][c] = data[c];
}

static void fill_y1_zedge(unimesh_patch_t* patch, real_t* data)
{
  DECLARE_UNIMESH_ZEDGE_ARRAY(a, patch);
  for (int ii = 0; ii <= patch->nx; ++ii)
    for (int kk = 0; kk < patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][0][kk][c] = data[c];
}

static void fill_y2_zedge(unimesh_patch_t* patch, real_t* data)
{
  DECLARE_UNIMESH_ZEDGE_ARRAY(a, patch);
  for (int ii = 0; ii <= patch->nx; ++ii)
    for (int kk = 0; kk < patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][patch->ny][kk][c] = data[c];
}

static void fill_z1_zedge(unimesh_patch_t* patch, real_t* data)
{
}

static void fill_z2_zedge(unimesh_patch_t* patch, real_t* data)
{
}

static void fill_x1_node(unimesh_patch_t* patch, real_t* data)
{
  DECLARE_UNIMESH_NODE_ARRAY(a, patch);
  for (int jj = 0; jj <= patch->ny; ++jj)
    for (int kk = 0; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[0][jj][kk][c] = data[c];
}

static void fill_x2_node(unimesh_patch_t* patch, real_t* data)
{
  DECLARE_UNIMESH_NODE_ARRAY(a, patch);
  for (int jj = 0; jj <= patch->ny; ++jj)
    for (int kk = 0; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[patch->nx][jj][kk][c] = data[c];
}

static void fill_y1_node(unimesh_patch_t* patch, real_t* data)
{
  DECLARE_UNIMESH_NODE_ARRAY(a, patch);
  for (int ii = 0; ii <= patch->nx; ++ii)
    for (int kk = 0; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][0][kk][c] = data[c];
}

static void fill_y2_node(unimesh_patch_t* patch, real_t* data)
{
  DECLARE_UNIMESH_NODE_ARRAY(a, patch);
  for (int ii = 0; ii <= patch->nx; ++ii)
    for (int kk = 0; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][patch->ny][kk][c] = data[c];
}

static void fill_z1_node(unimesh_patch_t* patch, real_t* data)
{
  DECLARE_UNIMESH_NODE_ARRAY(a, patch);
  for (int ii = 0; ii <= patch->nx; ++ii)
    for (int jj = 0; jj <= patch->ny; ++jj)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][jj][0][c] = data[c];
}

static void fill_z2_node(unimesh_patch_t* patch, real_t* data)
{
  DECLARE_UNIMESH_NODE_ARRAY(a, patch);
  for (int ii = 0; ii <= patch->nx; ++ii)
    for (int jj = 0; jj <= patch->ny; ++jj)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][jj][patch->nz][c] = data[c];
}

typedef void (*buffer_fill_func)(unimesh_patch_t* patch, real_t* data);

void unimesh_patch_fill_boundary(unimesh_patch_t* patch,
                                 unimesh_boundary_t boundary,
                                 real_t* data)
{
  static buffer_fill_func fill[8][6] =
  {
    {fill_x1_cell, fill_x2_cell, fill_y1_cell, fill_y2_cell, fill_z1_cell, fill_z2_cell},
    {fill_x1_xface, fill_x2_xface, fill_y1_xface, fill_y2_xface, fill_z1_xface, fill_z2_xface},
    {fill_x1_yface, fill_x2_yface, fill_y1_yface, fill_y2_yface, fill_z1_yface, fill_z2_yface},
    {fill_x1_zface, fill_x2_zface, fill_y1_zface, fill_y2_zface, fill_z1_zface, fill_z2_zface},
    {fill_x1_xedge, fill_x2_xedge, fill_y1_xedge, fill_y2_xedge, fill_z1_xedge, fill_z2_xedge},
    {fill_x1_yedge, fill_x2_yedge, fill_y1_yedge, fill_y2_yedge, fill_z1_yedge, fill_z2_yedge},
    {fill_x1_zedge, fill_x2_zedge, fill_y1_zedge, fill_y2_zedge, fill_z1_zedge, fill_z2_zedge},
    {fill_x1_node, fill_x2_node, fill_y1_node, fill_y2_node, fill_z1_node, fill_z2_node}
  };
  int c = (int)patch->centering;
  int b = (int)boundary;
  fill[c][b](patch, data);
}

