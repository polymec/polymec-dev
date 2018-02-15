// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

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
  size_t size = unimesh_patch_data_size(patch->centering, patch->nx, patch->ny, patch->nz, patch->nc);
  ASSERT(size == unimesh_patch_data_size(dest->centering, dest->nx, dest->ny, dest->nz, dest->nc));
  memcpy(dest->data, patch->data, size);
}

