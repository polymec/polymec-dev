// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/enumerable.h"
#include "geometry/unimesh_patch.h"

size_t unimesh_patch_data_size(unimesh_centering_t centering,
                               int nx, int ny, int nz, int nc)
{
  ASSERT(nx > 0);
  ASSERT(ny > 0);
  ASSERT(nz > 0);
  ASSERT(nc > 0);
  int num_data = 0;
  switch (centering)
  {
    case UNIMESH_NODE:  num_data = nc * (nx+1) * (ny+1) * (nz+1); break;
    case UNIMESH_XEDGE: num_data = nc * nx * (ny+1) * (nz+1); break;
    case UNIMESH_YEDGE: num_data = nc * (nx+1) * ny * (nz+1); break;
    case UNIMESH_ZEDGE: num_data = nc * (nx+1) * (ny+1) * nz; break;
    case UNIMESH_XFACE: num_data = nc * (nx+1) * ny * nz; break;
    case UNIMESH_YFACE: num_data = nc * nx * (ny+1) * nz; break;
    case UNIMESH_ZFACE: num_data = nc * nx * ny * (nz+1); break;
    case UNIMESH_CELL:  num_data = nc * (nx+2) * (ny+2) * (nz+2);
  }
  return sizeof(real_t) * num_data;
}

unimesh_patch_t* unimesh_patch_new(unimesh_centering_t centering,
                                   int nx, int ny, int nz, int nc)
{
  ASSERT(nx > 0);
  ASSERT(ny > 0);
  ASSERT(nz > 0);
  ASSERT(nc > 0);

  // We allocate one big slab of memory for storage and lean on C99's
  // VLA semantics.
  size_t data_size = unimesh_patch_data_size(centering, nx, ny, nz, nc);
  size_t storage_size = sizeof(unimesh_patch_t) + data_size;
  unimesh_patch_t* p = polymec_malloc(storage_size);
  p->data = (char*)p + sizeof(unimesh_patch_t);
  memset(p->data, 0, data_size);
  p->nx = nx;
  p->ny = ny;
  p->nz = nz;
  p->nc = nc;
  p->centering = centering;
  return p;
}

unimesh_patch_t* unimesh_patch_with_buffer(unimesh_centering_t centering,
                                           int nx, int ny, int nz, int nc,
                                           void* buffer)
{
  ASSERT(nx > 0);
  ASSERT(ny > 0);
  ASSERT(nz > 0);
  ASSERT(nc > 0);

  unimesh_patch_t* p = polymec_malloc(sizeof(unimesh_patch_t));
  p->data = buffer;
  p->nx = nx;
  p->ny = ny;
  p->nz = nz;
  p->nc = nc;
  p->centering = centering;
  return p;
}

unimesh_patch_t* unimesh_patch_clone(unimesh_patch_t* patch)
{
  unimesh_patch_t* clone = unimesh_patch_new(patch->centering,
                                             patch->nx, patch->ny, patch->nz,
                                             patch->nc);
  unimesh_patch_copy(patch, clone);
  return clone;
}

void unimesh_patch_free(unimesh_patch_t* patch)
{
  polymec_free(patch);
}

void unimesh_patch_copy(unimesh_patch_t* patch,
                        unimesh_patch_t* dest)
{
  size_t size = unimesh_patch_data_size(patch->centering, patch->nx, patch->ny,
                                        patch->nz, patch->nc);
  ASSERT(size == unimesh_patch_data_size(dest->centering, dest->nx, dest->ny,
                                         dest->nz, dest->nc));
  memcpy(dest->data, patch->data, size);
}

void unimesh_patch_get_box(unimesh_patch_t* patch,
                           unimesh_patch_box_t* box)
{
  if (patch->centering == UNIMESH_CELL)
  {
    box->i1 = 1;
    box->i2 = patch->nx+1;
    box->j1 = 1;
    box->j2 = patch->ny+1;
    box->k1 = 1;
    box->k2 = patch->nz+1;
  }
  else
  {
    box->i1 = 0;
    box->j1 = 0;
    box->k1 = 0;
    if (patch->centering == UNIMESH_XFACE)
    {
      box->i2 = patch->nx+1;
      box->j2 = patch->ny;
      box->k2 = patch->nz;
    }
    else if (patch->centering == UNIMESH_YFACE)
    {
      box->i2 = patch->nx;
      box->j2 = patch->ny+1;
      box->k2 = patch->nz;
    }
    else if (patch->centering == UNIMESH_ZFACE)
    {
      box->i2 = patch->nx;
      box->j2 = patch->ny;
      box->k2 = patch->nz+1;
    }
    else if (patch->centering == UNIMESH_XEDGE)
    {
      box->i2 = patch->nx;
      box->j2 = patch->ny+1;
      box->k2 = patch->nz+1;
    }
    else if (patch->centering == UNIMESH_YEDGE)
    {
      box->i2 = patch->nx+1;
      box->j2 = patch->ny;
      box->k2 = patch->nz+1;
    }
    else if (patch->centering == UNIMESH_ZEDGE)
    {
      box->i2 = patch->nx+1;
      box->j2 = patch->ny+1;
      box->k2 = patch->nz;
    }
    else if (patch->centering == UNIMESH_NODE)
    {
      box->i2 = patch->nx+1;
      box->j2 = patch->ny+1;
      box->k2 = patch->nz+1;
    }
  }
}

void unimesh_patch_get_boundary_box(unimesh_patch_t* patch,
                                    unimesh_boundary_t boundary,
                                    unimesh_patch_box_t* box)
{
  unimesh_patch_get_box(patch, box);

  if (patch->centering == UNIMESH_CELL)
  {
    switch (boundary)
    {
      case UNIMESH_X1_BOUNDARY:
        box->i2 = box->i1 + 1;
        unimesh_patch_box_shift(box, -1, 0, 0);
        break;
      case UNIMESH_X2_BOUNDARY:
        box->i1 = box->i2 - 1;
        unimesh_patch_box_shift(box, 1, 0, 0);
        break;
      case UNIMESH_Y1_BOUNDARY:
        box->j2 = box->j1 + 1;
        unimesh_patch_box_shift(box, 0, -1, 0);
        break;
      case UNIMESH_Y2_BOUNDARY:
        box->j1 = box->j2 - 1;
        unimesh_patch_box_shift(box, 0, 1, 0);
        break;
      case UNIMESH_Z1_BOUNDARY:
        box->k2 = box->k1 + 1;
        unimesh_patch_box_shift(box, 0, 0, -1);
        break;
      case UNIMESH_Z2_BOUNDARY:
        box->k1 = box->k2 - 1;
        unimesh_patch_box_shift(box, 0, 0, 1);
    }
  }
  else
  {
    switch (boundary)
    {
      case UNIMESH_X1_BOUNDARY:
        box->i2 = box->i1 + 1;
        break;
      case UNIMESH_X2_BOUNDARY:
        box->i1 = box->i2 - 1;
        break;
      case UNIMESH_Y1_BOUNDARY:
        box->j2 = box->j1 + 1;
        break;
      case UNIMESH_Y2_BOUNDARY:
        box->j1 = box->j2 - 1;
        break;
      case UNIMESH_Z1_BOUNDARY:
        box->k2 = box->k1 + 1;
        break;
      case UNIMESH_Z2_BOUNDARY:
        box->k1 = box->k2 - 1;
    }
  }
}

void unimesh_patch_box_bisect(unimesh_patch_box_t* box, int axis, int half)
{
  if (axis == 0)
  {
    if (half == 0)
      box->i2 /= 2;
    else
      box->i1 = box->i2/2;
  }
  else if (axis == 1)
  {
    if (half == 0)
      box->j2 /= 2;
    else
      box->j1 = box->j2/2;
  }
  else
  {
    ASSERT(axis == 2);
    if (half == 0)
      box->k2 /= 2;
    else
      box->k1 = box->k2/2;
  }
}

real_enumerable_generator_t* unimesh_patch_enumerate(unimesh_patch_t* patch)
{
  size_t num_values = unimesh_patch_data_size(patch->centering, patch->nx, patch->ny, patch->nz, patch->nc) / sizeof(real_t);
  return real_enumerable_generator_from_array((real_t*)patch->data, num_values, NULL);
}

