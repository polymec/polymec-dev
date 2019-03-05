// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "geometry/unimesh.h"
#include "geometry/unimesh_patch.h"

static void copy_x1_cell_to(unimesh_patch_t* patch, void* buffer)
{
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->ny+2, patch->nz+2, patch->nc);
  DECLARE_UNIMESH_CELL_ARRAY(a, patch);
  for (int jj = 1; jj <= patch->ny; ++jj)
    for (int kk = 1; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        buf[jj][kk][c] = a[1][jj][kk][c];
}

static void copy_x2_cell_to(unimesh_patch_t* patch, void* buffer)
{
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->ny+2, patch->nz+2, patch->nc);
  DECLARE_UNIMESH_CELL_ARRAY(a, patch);
  for (int jj = 1; jj <= patch->ny; ++jj)
    for (int kk = 1; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        buf[jj][kk][c] = a[patch->nx][jj][kk][c];
}

static void copy_y1_cell_to(unimesh_patch_t* patch, void* buffer)
{
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->nx+2, patch->nz+2, patch->nc);
  DECLARE_UNIMESH_CELL_ARRAY(a, patch);
  for (int ii = 1; ii <= patch->nx; ++ii)
    for (int kk = 1; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        buf[ii][kk][c] = a[ii][1][kk][c];
}

static void copy_y2_cell_to(unimesh_patch_t* patch, void* buffer)
{
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->nx+2, patch->nz+2, patch->nc);
  DECLARE_UNIMESH_CELL_ARRAY(a, patch);
  for (int ii = 1; ii <= patch->nx; ++ii)
    for (int kk = 1; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        buf[ii][kk][c] = a[ii][patch->ny][kk][c];
}

static void copy_z1_cell_to(unimesh_patch_t* patch, void* buffer)
{
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->nx+2, patch->ny+2, patch->nc);
  DECLARE_UNIMESH_CELL_ARRAY(a, patch);
  for (int ii = 1; ii <= patch->nx; ++ii)
    for (int jj = 1; jj <= patch->ny; ++jj)
      for (int c = 0; c < patch->nc; ++c)
        buf[ii][jj][c] = a[ii][jj][1][c];
}

static void copy_z2_cell_to(unimesh_patch_t* patch, void* buffer)
{
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->nx+2, patch->ny+2, patch->nc);
  DECLARE_UNIMESH_CELL_ARRAY(a, patch);
  for (int ii = 1; ii <= patch->nx; ++ii)
    for (int jj = 1; jj <= patch->ny; ++jj)
      for (int c = 0; c < patch->nc; ++c)
        buf[ii][jj][c] = a[ii][jj][patch->nz][c];
}

static void copy_x1_xface_to(unimesh_patch_t* patch, void* buffer)
{
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->ny, patch->nz, patch->nc);
  DECLARE_UNIMESH_XFACE_ARRAY(a, patch);
  for (int jj = 0; jj < patch->ny; ++jj)
    for (int kk = 0; kk < patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        buf[jj][kk][c] = a[0][jj][kk][c];
}

static void copy_x2_xface_to(unimesh_patch_t* patch, void* buffer)
{
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->ny, patch->nz, patch->nc);
  DECLARE_UNIMESH_XFACE_ARRAY(a, patch);
  for (int jj = 0; jj < patch->ny; ++jj)
    for (int kk = 0; kk < patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        buf[jj][kk][c] = a[patch->nx][jj][kk][c];
}

static void copy_y1_xface_to(unimesh_patch_t* patch, void* buffer)
{
}

static void copy_y2_xface_to(unimesh_patch_t* patch, void* buffer)
{
}

static void copy_z1_xface_to(unimesh_patch_t* patch, void* buffer)
{
}

static void copy_z2_xface_to(unimesh_patch_t* patch, void* buffer)
{
}

static void copy_x1_yface_to(unimesh_patch_t* patch, void* buffer)
{
}

static void copy_x2_yface_to(unimesh_patch_t* patch, void* buffer)
{
}

static void copy_y1_yface_to(unimesh_patch_t* patch, void* buffer)
{
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->nx, patch->nz, patch->nc);
  DECLARE_UNIMESH_YFACE_ARRAY(a, patch);
  for (int ii = 0; ii < patch->nx; ++ii)
    for (int kk = 0; kk < patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        buf[ii][kk][c] = a[ii][0][kk][c];
}

static void copy_y2_yface_to(unimesh_patch_t* patch, void* buffer)
{
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->nx, patch->nz, patch->nc);
  DECLARE_UNIMESH_YFACE_ARRAY(a, patch);
  for (int ii = 0; ii < patch->nx; ++ii)
    for (int kk = 0; kk < patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        buf[ii][kk][c] = a[ii][patch->ny][kk][c];
}

static void copy_z1_yface_to(unimesh_patch_t* patch, void* buffer)
{
}

static void copy_z2_yface_to(unimesh_patch_t* patch, void* buffer)
{
}

static void copy_x1_zface_to(unimesh_patch_t* patch, void* buffer)
{
}

static void copy_x2_zface_to(unimesh_patch_t* patch, void* buffer)
{
}

static void copy_y1_zface_to(unimesh_patch_t* patch, void* buffer)
{
}

static void copy_y2_zface_to(unimesh_patch_t* patch, void* buffer)
{
}

static void copy_z1_zface_to(unimesh_patch_t* patch, void* buffer)
{
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->nx, patch->ny, patch->nc);
  DECLARE_UNIMESH_ZFACE_ARRAY(a, patch);
  for (int ii = 0; ii < patch->nx; ++ii)
    for (int jj = 0; jj < patch->ny; ++jj)
      for (int c = 0; c < patch->nc; ++c)
        buf[ii][jj][c] = a[ii][jj][0][c];
}

static void copy_z2_zface_to(unimesh_patch_t* patch, void* buffer)
{
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->nx, patch->ny, patch->nc);
  DECLARE_UNIMESH_ZFACE_ARRAY(a, patch);
  for (int ii = 0; ii < patch->nx; ++ii)
    for (int jj = 0; jj < patch->ny; ++jj)
      for (int c = 0; c < patch->nc; ++c)
        buf[ii][jj][c] = a[ii][jj][patch->nz][c];
}

static void copy_x1_xedge_to(unimesh_patch_t* patch, void* buffer)
{
}

static void copy_x2_xedge_to(unimesh_patch_t* patch, void* buffer)
{
}

static void copy_y1_xedge_to(unimesh_patch_t* patch, void* buffer)
{
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->nx, patch->nz+1, patch->nc);
  DECLARE_UNIMESH_XEDGE_ARRAY(a, patch);
  for (int ii = 0; ii < patch->nx; ++ii)
    for (int kk = 0; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        buf[ii][kk][c] = a[ii][0][kk][c];
}

static void copy_y2_xedge_to(unimesh_patch_t* patch, void* buffer)
{
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->nx, patch->nz+1, patch->nc);
  DECLARE_UNIMESH_XEDGE_ARRAY(a, patch);
  for (int ii = 0; ii < patch->nx; ++ii)
    for (int kk = 0; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        buf[ii][kk][c] = a[ii][patch->ny][kk][c];
}

static void copy_z1_xedge_to(unimesh_patch_t* patch, void* buffer)
{
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->nx, patch->ny+1, patch->nc);
  DECLARE_UNIMESH_XEDGE_ARRAY(a, patch);
  for (int ii = 0; ii < patch->nx; ++ii)
    for (int jj = 0; jj <= patch->ny; ++jj)
      for (int c = 0; c < patch->nc; ++c)
        buf[ii][jj][c] = a[ii][jj][0][c];
}

static void copy_z2_xedge_to(unimesh_patch_t* patch, void* buffer)
{
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->nx, patch->ny+1, patch->nc);
  DECLARE_UNIMESH_XEDGE_ARRAY(a, patch);
  for (int ii = 0; ii < patch->nx; ++ii)
    for (int jj = 0; jj <= patch->ny; ++jj)
      for (int c = 0; c < patch->nc; ++c)
        buf[ii][jj][c] = a[ii][jj][patch->nz][c];
}

static void copy_x1_yedge_to(unimesh_patch_t* patch, void* buffer)
{
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->ny, patch->nz+1, patch->nc);
  DECLARE_UNIMESH_YEDGE_ARRAY(a, patch);
  for (int jj = 0; jj < patch->ny; ++jj)
    for (int kk = 0; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        buf[jj][kk][c] = a[0][jj][kk][c];
}

static void copy_x2_yedge_to(unimesh_patch_t* patch, void* buffer)
{
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->ny, patch->nz+1, patch->nc);
  DECLARE_UNIMESH_YEDGE_ARRAY(a, patch);
  for (int jj = 0; jj < patch->ny; ++jj)
    for (int kk = 0; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        buf[jj][kk][c] = a[patch->nx][jj][kk][c];
}

static void copy_y1_yedge_to(unimesh_patch_t* patch, void* buffer)
{
}

static void copy_y2_yedge_to(unimesh_patch_t* patch, void* buffer)
{
}

static void copy_z1_yedge_to(unimesh_patch_t* patch, void* buffer)
{
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->nx+1, patch->ny, patch->nc);
  DECLARE_UNIMESH_YEDGE_ARRAY(a, patch);
  for (int ii = 0; ii <= patch->nx; ++ii)
    for (int jj = 0; jj < patch->ny; ++jj)
      for (int c = 0; c < patch->nc; ++c)
        buf[ii][jj][c] = a[ii][jj][0][c];
}

static void copy_z2_yedge_to(unimesh_patch_t* patch, void* buffer)
{
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->nx+1, patch->ny, patch->nc);
  DECLARE_UNIMESH_YEDGE_ARRAY(a, patch);
  for (int ii = 0; ii <= patch->nx; ++ii)
    for (int jj = 0; jj < patch->ny; ++jj)
      for (int c = 0; c < patch->nc; ++c)
        buf[ii][jj][c] = a[ii][jj][patch->nz][c];
}

static void copy_x1_zedge_to(unimesh_patch_t* patch, void* buffer)
{
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->ny+1, patch->nz, patch->nc);
  DECLARE_UNIMESH_ZEDGE_ARRAY(a, patch);
  for (int jj = 0; jj <= patch->ny; ++jj)
    for (int kk = 0; kk < patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        buf[jj][kk][c] = a[0][jj][kk][c];
}

static void copy_x2_zedge_to(unimesh_patch_t* patch, void* buffer)
{
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->ny+1, patch->nz, patch->nc);
  DECLARE_UNIMESH_ZEDGE_ARRAY(a, patch);
  for (int jj = 0; jj <= patch->ny; ++jj)
    for (int kk = 0; kk < patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        buf[jj][kk][c] = a[patch->nx][jj][kk][c];
}

static void copy_y1_zedge_to(unimesh_patch_t* patch, void* buffer)
{
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->nx+1, patch->nz, patch->nc);
  DECLARE_UNIMESH_ZEDGE_ARRAY(a, patch);
  for (int ii = 0; ii <= patch->nx; ++ii)
    for (int kk = 0; kk < patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        buf[ii][kk][c] = a[ii][0][kk][c];
}

static void copy_y2_zedge_to(unimesh_patch_t* patch, void* buffer)
{
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->nx+1, patch->nz, patch->nc);
  DECLARE_UNIMESH_ZEDGE_ARRAY(a, patch);
  for (int ii = 0; ii <= patch->nx; ++ii)
    for (int kk = 0; kk < patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        buf[ii][kk][c] = a[ii][patch->ny][kk][c];
}

static void copy_z1_zedge_to(unimesh_patch_t* patch, void* buffer)
{
}

static void copy_z2_zedge_to(unimesh_patch_t* patch, void* buffer)
{
}

static void copy_x1_node_to(unimesh_patch_t* patch, void* buffer)
{
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->ny+1, patch->nz+1, patch->nc);
  DECLARE_UNIMESH_NODE_ARRAY(a, patch);
  for (int jj = 0; jj <= patch->ny; ++jj)
    for (int kk = 0; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        buf[jj][kk][c] = a[0][jj][kk][c];
}

static void copy_x2_node_to(unimesh_patch_t* patch, void* buffer)
{
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->ny+1, patch->nz+1, patch->nc);
  DECLARE_UNIMESH_NODE_ARRAY(a, patch);
  for (int jj = 0; jj <= patch->ny; ++jj)
    for (int kk = 0; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        buf[jj][kk][c] = a[patch->nx][jj][kk][c];
}

static void copy_y1_node_to(unimesh_patch_t* patch, void* buffer)
{
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->nx+1, patch->nz+1, patch->nc);
  DECLARE_UNIMESH_NODE_ARRAY(a, patch);
  for (int ii = 0; ii <= patch->nx; ++ii)
    for (int kk = 0; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        buf[ii][kk][c] = a[ii][0][kk][c];
}

static void copy_y2_node_to(unimesh_patch_t* patch, void* buffer)
{
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->nx+1, patch->nz+1, patch->nc);
  DECLARE_UNIMESH_NODE_ARRAY(a, patch);
  for (int ii = 0; ii <= patch->nx; ++ii)
    for (int kk = 0; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        buf[ii][kk][c] = a[ii][patch->ny][kk][c];
}

static void copy_z1_node_to(unimesh_patch_t* patch, void* buffer)
{
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->nx+1, patch->ny+1, patch->nc);
  DECLARE_UNIMESH_NODE_ARRAY(a, patch);
  for (int ii = 0; ii <= patch->nx; ++ii)
    for (int jj = 0; jj <= patch->ny; ++jj)
      for (int c = 0; c < patch->nc; ++c)
        buf[ii][jj][c] = a[ii][jj][0][c];
}

static void copy_z2_node_to(unimesh_patch_t* patch, void* buffer)
{
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->nx+1, patch->ny+1, patch->nc);
  DECLARE_UNIMESH_NODE_ARRAY(a, patch);
  for (int ii = 0; ii <= patch->nx; ++ii)
    for (int jj = 0; jj <= patch->ny; ++jj)
      for (int c = 0; c < patch->nc; ++c)
        buf[ii][jj][c] = a[ii][jj][patch->nz][c];
}

static void copy_x1_cell_from(unimesh_patch_t* patch, void* buffer)
{
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->ny+2, patch->nz+2, patch->nc);
  DECLARE_UNIMESH_CELL_ARRAY(a, patch);
  for (int jj = 1; jj <= patch->ny; ++jj)
    for (int kk = 1; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[0][jj][kk][c] = buf[jj][kk][c];
}

static void copy_x2_cell_from(unimesh_patch_t* patch, void* buffer)
{
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->ny+2, patch->nz+2, patch->nc);
  DECLARE_UNIMESH_CELL_ARRAY(a, patch);
  for (int jj = 1; jj <= patch->ny; ++jj)
    for (int kk = 1; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[patch->nx+1][jj][kk][c] = buf[jj][kk][c];
}

static void copy_y1_cell_from(unimesh_patch_t* patch, void* buffer)
{
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->nx+2, patch->nz+2, patch->nc);
  DECLARE_UNIMESH_CELL_ARRAY(a, patch);
  for (int ii = 1; ii <= patch->nx; ++ii)
    for (int kk = 1; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][0][kk][c] = buf[ii][kk][c];
}

static void copy_y2_cell_from(unimesh_patch_t* patch, void* buffer)
{
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->nx+2, patch->nz+2, patch->nc);
  DECLARE_UNIMESH_CELL_ARRAY(a, patch);
  for (int ii = 1; ii <= patch->nx; ++ii)
    for (int kk = 1; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][patch->ny+1][kk][c] = buf[ii][kk][c];
}

static void copy_z1_cell_from(unimesh_patch_t* patch, void* buffer)
{
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->nx+2, patch->ny+2, patch->nc);
  DECLARE_UNIMESH_CELL_ARRAY(a, patch);
  for (int ii = 1; ii <= patch->nx; ++ii)
    for (int jj = 1; jj <= patch->ny; ++jj)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][jj][0][c] = buf[ii][jj][c];
}

static void copy_z2_cell_from(unimesh_patch_t* patch, void* buffer)
{
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->nx+2, patch->ny+2, patch->nc);
  DECLARE_UNIMESH_CELL_ARRAY(a, patch);
  for (int ii = 1; ii <= patch->nx; ++ii)
    for (int jj = 1; jj <= patch->ny; ++jj)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][jj][patch->nz+1][c] = buf[ii][jj][c];
}

static void copy_x1_xface_from(unimesh_patch_t* patch, void* buffer)
{
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->ny, patch->nz, patch->nc);
  DECLARE_UNIMESH_XFACE_ARRAY(a, patch);
  for (int jj = 0; jj < patch->ny; ++jj)
    for (int kk = 0; kk < patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[0][jj][kk][c] = buf[jj][kk][c];
}

static void copy_x2_xface_from(unimesh_patch_t* patch, void* buffer)
{
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->ny, patch->nz, patch->nc);
  DECLARE_UNIMESH_XFACE_ARRAY(a, patch);
  for (int jj = 0; jj < patch->ny; ++jj)
    for (int kk = 0; kk < patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[patch->nx][jj][kk][c] = buf[jj][kk][c];
}

static void copy_y1_xface_from(unimesh_patch_t* patch, void* buffer)
{
}

static void copy_y2_xface_from(unimesh_patch_t* patch, void* buffer)
{
}

static void copy_z1_xface_from(unimesh_patch_t* patch, void* buffer)
{
}

static void copy_z2_xface_from(unimesh_patch_t* patch, void* buffer)
{
}

static void copy_x1_yface_from(unimesh_patch_t* patch, void* buffer)
{
}

static void copy_x2_yface_from(unimesh_patch_t* patch, void* buffer)
{
}

static void copy_y1_yface_from(unimesh_patch_t* patch, void* buffer)
{
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->nx, patch->nz, patch->nc);
  DECLARE_UNIMESH_YFACE_ARRAY(a, patch);
  for (int ii = 0; ii < patch->nx; ++ii)
    for (int kk = 0; kk < patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][0][kk][c] = buf[ii][kk][c];
}

static void copy_y2_yface_from(unimesh_patch_t* patch, void* buffer)
{
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->nx, patch->nz, patch->nc);
  DECLARE_UNIMESH_YFACE_ARRAY(a, patch);
  for (int ii = 0; ii < patch->nx; ++ii)
    for (int kk = 0; kk < patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][patch->ny][kk][c] = buf[ii][kk][c];
}

static void copy_z1_yface_from(unimesh_patch_t* patch, void* buffer)
{
}

static void copy_z2_yface_from(unimesh_patch_t* patch, void* buffer)
{
}

static void copy_x1_zface_from(unimesh_patch_t* patch, void* buffer)
{
}

static void copy_x2_zface_from(unimesh_patch_t* patch, void* buffer)
{
}

static void copy_y1_zface_from(unimesh_patch_t* patch, void* buffer)
{
}

static void copy_y2_zface_from(unimesh_patch_t* patch, void* buffer)
{
}

static void copy_z1_zface_from(unimesh_patch_t* patch, void* buffer)
{
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->nx, patch->ny, patch->nc);
  DECLARE_UNIMESH_ZFACE_ARRAY(a, patch);
  for (int ii = 0; ii < patch->nx; ++ii)
    for (int jj = 0; jj < patch->ny; ++jj)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][jj][0][c] = buf[ii][jj][c];
}

static void copy_z2_zface_from(unimesh_patch_t* patch, void* buffer)
{
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->nx, patch->ny, patch->nc);
  DECLARE_UNIMESH_ZFACE_ARRAY(a, patch);
  for (int ii = 0; ii < patch->nx; ++ii)
    for (int jj = 0; jj < patch->ny; ++jj)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][jj][patch->nz][c] = buf[ii][jj][c];
}

static void copy_x1_xedge_from(unimesh_patch_t* patch, void* buffer)
{
}

static void copy_x2_xedge_from(unimesh_patch_t* patch, void* buffer)
{
}

static void copy_y1_xedge_from(unimesh_patch_t* patch, void* buffer)
{
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->nx, patch->nz+1, patch->nc);
  DECLARE_UNIMESH_XEDGE_ARRAY(a, patch);
  for (int ii = 0; ii < patch->nx; ++ii)
    for (int kk = 0; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][0][kk][c] = buf[ii][kk][c];
}

static void copy_y2_xedge_from(unimesh_patch_t* patch, void* buffer)
{
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->nx, patch->nz+1, patch->nc);
  DECLARE_UNIMESH_XEDGE_ARRAY(a, patch);
  for (int ii = 0; ii < patch->nx; ++ii)
    for (int kk = 0; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][patch->ny][kk][c] = buf[ii][kk][c];
}

static void copy_z1_xedge_from(unimesh_patch_t* patch, void* buffer)
{
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->nx, patch->ny+1, patch->nc);
  DECLARE_UNIMESH_XEDGE_ARRAY(a, patch);
  for (int ii = 0; ii < patch->nx; ++ii)
    for (int jj = 0; jj <= patch->ny; ++jj)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][jj][0][c] = buf[ii][jj][c];
}

static void copy_z2_xedge_from(unimesh_patch_t* patch, void* buffer)
{
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->nx, patch->ny+1, patch->nc);
  DECLARE_UNIMESH_XEDGE_ARRAY(a, patch);
  for (int ii = 0; ii < patch->nx; ++ii)
    for (int jj = 0; jj <= patch->ny; ++jj)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][jj][patch->nz][c] = buf[ii][jj][c];
}

static void copy_x1_yedge_from(unimesh_patch_t* patch, void* buffer)
{
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->ny, patch->nz+1, patch->nc);
  DECLARE_UNIMESH_YEDGE_ARRAY(a, patch);
  for (int jj = 0; jj < patch->ny; ++jj)
    for (int kk = 0; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[0][jj][kk][c] = buf[jj][kk][c];
}

static void copy_x2_yedge_from(unimesh_patch_t* patch, void* buffer)
{
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->ny, patch->nz+1, patch->nc);
  DECLARE_UNIMESH_YEDGE_ARRAY(a, patch);
  for (int jj = 0; jj < patch->ny; ++jj)
    for (int kk = 0; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[patch->nx][jj][kk][c] = buf[jj][kk][c];
}

static void copy_y1_yedge_from(unimesh_patch_t* patch, void* buffer)
{
}

static void copy_y2_yedge_from(unimesh_patch_t* patch, void* buffer)
{
}

static void copy_z1_yedge_from(unimesh_patch_t* patch, void* buffer)
{
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->nx+1, patch->ny, patch->nc);
  DECLARE_UNIMESH_YEDGE_ARRAY(a, patch);
  for (int ii = 0; ii <= patch->nx; ++ii)
    for (int jj = 0; jj < patch->ny; ++jj)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][jj][0][c] = buf[ii][jj][c];
}

static void copy_z2_yedge_from(unimesh_patch_t* patch, void* buffer)
{
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->nx+1, patch->ny, patch->nc);
  DECLARE_UNIMESH_YEDGE_ARRAY(a, patch);
  for (int ii = 0; ii <= patch->nx; ++ii)
    for (int jj = 0; jj < patch->ny; ++jj)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][jj][patch->nz][c] = buf[ii][jj][c];
}

static void copy_x1_zedge_from(unimesh_patch_t* patch, void* buffer)
{
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->ny+1, patch->nz, patch->nc);
  DECLARE_UNIMESH_ZEDGE_ARRAY(a, patch);
  for (int jj = 0; jj <= patch->ny; ++jj)
    for (int kk = 0; kk < patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[0][jj][kk][c] = buf[jj][kk][c];
}

static void copy_x2_zedge_from(unimesh_patch_t* patch, void* buffer)
{
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->ny+1, patch->nz, patch->nc);
  DECLARE_UNIMESH_ZEDGE_ARRAY(a, patch);
  for (int jj = 0; jj <= patch->ny; ++jj)
    for (int kk = 0; kk < patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[patch->nx][jj][kk][c] = buf[jj][kk][c];
}

static void copy_y1_zedge_from(unimesh_patch_t* patch, void* buffer)
{
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->nx+1, patch->nz, patch->nc);
  DECLARE_UNIMESH_ZEDGE_ARRAY(a, patch);
  for (int ii = 0; ii <= patch->nx; ++ii)
    for (int kk = 0; kk < patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][0][kk][c] = buf[ii][kk][c];
}

static void copy_y2_zedge_from(unimesh_patch_t* patch, void* buffer)
{
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->nx+1, patch->nz, patch->nc);
  DECLARE_UNIMESH_ZEDGE_ARRAY(a, patch);
  for (int ii = 0; ii <= patch->nx; ++ii)
    for (int kk = 0; kk < patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][patch->ny][kk][c] = buf[ii][kk][c];
}

static void copy_z1_zedge_from(unimesh_patch_t* patch, void* buffer)
{
}

static void copy_z2_zedge_from(unimesh_patch_t* patch, void* buffer)
{
}

static void copy_x1_node_from(unimesh_patch_t* patch, void* buffer)
{
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->ny+1, patch->nz+1, patch->nc);
  DECLARE_UNIMESH_NODE_ARRAY(a, patch);
  for (int jj = 0; jj <= patch->ny; ++jj)
    for (int kk = 0; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[0][jj][kk][c] = buf[jj][kk][c];
}

static void copy_x2_node_from(unimesh_patch_t* patch, void* buffer)
{
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->ny+1, patch->nz+1, patch->nc);
  DECLARE_UNIMESH_NODE_ARRAY(a, patch);
  for (int jj = 0; jj <= patch->ny; ++jj)
    for (int kk = 0; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[patch->nx][jj][kk][c] = buf[jj][kk][c];
}

static void copy_y1_node_from(unimesh_patch_t* patch, void* buffer)
{
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->nx+1, patch->nz+1, patch->nc);
  DECLARE_UNIMESH_NODE_ARRAY(a, patch);
  for (int ii = 0; ii <= patch->nx; ++ii)
    for (int kk = 0; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][0][kk][c] = buf[ii][kk][c];
}

static void copy_y2_node_from(unimesh_patch_t* patch, void* buffer)
{
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->nx+1, patch->nz+1, patch->nc);
  DECLARE_UNIMESH_NODE_ARRAY(a, patch);
  for (int ii = 0; ii <= patch->nx; ++ii)
    for (int kk = 0; kk <= patch->nz; ++kk)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][patch->ny][kk][c] = buf[ii][kk][c];
}

static void copy_z1_node_from(unimesh_patch_t* patch, void* buffer)
{
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->nx+1, patch->ny+1, patch->nc);
  DECLARE_UNIMESH_NODE_ARRAY(a, patch);
  for (int ii = 0; ii <= patch->nx; ++ii)
    for (int jj = 0; jj <= patch->ny; ++jj)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][jj][0][c] = buf[ii][jj][c];
}

static void copy_z2_node_from(unimesh_patch_t* patch, void* buffer)
{
  DECLARE_3D_ARRAY(real_t, buf, buffer, patch->nx+1, patch->ny+1, patch->nc);
  DECLARE_UNIMESH_NODE_ARRAY(a, patch);
  for (int ii = 0; ii <= patch->nx; ++ii)
    for (int jj = 0; jj <= patch->ny; ++jj)
      for (int c = 0; c < patch->nc; ++c)
        a[ii][jj][patch->nz][c] = buf[ii][jj][c];
}

typedef void (*buffer_copy_func)(unimesh_patch_t* patch, void* buffer);

void unimesh_patch_copy_bvalues_to_buffer(unimesh_patch_t* patch,
                                          unimesh_boundary_t boundary,
                                          void* buffer);
void unimesh_patch_copy_bvalues_to_buffer(unimesh_patch_t* patch,
                                          unimesh_boundary_t boundary,
                                          void* buffer)
{
  static buffer_copy_func copy_to[8][6] =
  {
    {copy_x1_cell_to, copy_x2_cell_to, copy_y1_cell_to, copy_y2_cell_to, copy_z1_cell_to, copy_z2_cell_to},
    {copy_x1_xface_to, copy_x2_xface_to, copy_y1_xface_to, copy_y2_xface_to, copy_z1_xface_to, copy_z2_xface_to},
    {copy_x1_yface_to, copy_x2_yface_to, copy_y1_yface_to, copy_y2_yface_to, copy_z1_yface_to, copy_z2_yface_to},
    {copy_x1_zface_to, copy_x2_zface_to, copy_y1_zface_to, copy_y2_zface_to, copy_z1_zface_to, copy_z2_zface_to},
    {copy_x1_xedge_to, copy_x2_xedge_to, copy_y1_xedge_to, copy_y2_xedge_to, copy_z1_xedge_to, copy_z2_xedge_to},
    {copy_x1_yedge_to, copy_x2_yedge_to, copy_y1_yedge_to, copy_y2_yedge_to, copy_z1_yedge_to, copy_z2_yedge_to},
    {copy_x1_zedge_to, copy_x2_zedge_to, copy_y1_zedge_to, copy_y2_zedge_to, copy_z1_zedge_to, copy_z2_zedge_to},
    {copy_x1_node_to, copy_x2_node_to, copy_y1_node_to, copy_y2_node_to, copy_z1_node_to, copy_z2_node_to}
  };
  int c = (int)patch->centering;
  int b = (int)boundary;
  copy_to[c][b](patch, buffer);
}

void unimesh_patch_copy_bvalues_from_buffer(unimesh_patch_t* patch,
                                            unimesh_boundary_t boundary,
                                            void* buffer);
void unimesh_patch_copy_bvalues_from_buffer(unimesh_patch_t* patch,
                                            unimesh_boundary_t boundary,
                                            void* buffer)
{
  static buffer_copy_func copy_from[8][6] =
  {
    {copy_x1_cell_from, copy_x2_cell_from, copy_y1_cell_from, copy_y2_cell_from, copy_z1_cell_from, copy_z2_cell_from},
    {copy_x1_xface_from, copy_x2_xface_from, copy_y1_xface_from, copy_y2_xface_from, copy_z1_xface_from, copy_z2_xface_from},
    {copy_x1_yface_from, copy_x2_yface_from, copy_y1_yface_from, copy_y2_yface_from, copy_z1_yface_from, copy_z2_yface_from},
    {copy_x1_zface_from, copy_x2_zface_from, copy_y1_zface_from, copy_y2_zface_from, copy_z1_zface_from, copy_z2_zface_from},
    {copy_x1_xedge_from, copy_x2_xedge_from, copy_y1_xedge_from, copy_y2_xedge_from, copy_z1_xedge_from, copy_z2_xedge_from},
    {copy_x1_yedge_from, copy_x2_yedge_from, copy_y1_yedge_from, copy_y2_yedge_from, copy_z1_yedge_from, copy_z2_yedge_from},
    {copy_x1_zedge_from, copy_x2_zedge_from, copy_y1_zedge_from, copy_y2_zedge_from, copy_z1_zedge_from, copy_z2_zedge_from},
    {copy_x1_node_from, copy_x2_node_from, copy_y1_node_from, copy_y2_node_from, copy_z1_node_from, copy_z2_node_from}
  };
  int c = (int)patch->centering;
  int b = (int)boundary;
  copy_from[c][b](patch, buffer);
}

