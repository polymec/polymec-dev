// Copyright (c) 2012-2014, Jeffrey N. Johnson
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this 
// list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice, 
// this list of conditions and the following disclaimer in the documentation 
// and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include <gc/gc.h>
#include "geometry/patch_grid.h"

static void patch_grid_free(void* ctx, void* dummy)
{
  patch_grid_t* grid = ctx;
  ptr_slist_free(grid->subgrids);
}

patch_grid_t* patch_grid_new(bbox_t* domain, int nx, int ny, int nz)
{
  ASSERT(domain->x1 < domain->x2);
  ASSERT(domain->y1 < domain->y2);
  ASSERT(domain->z1 < domain->z2);
  ASSERT(nx > 0);
  ASSERT(ny > 0);
  ASSERT(nz > 0);
  patch_grid_t* grid = GC_MALLOC(sizeof(patch_grid_t));
  grid->domain = *domain;
  grid->nx = nx;
  grid->ny = ny;
  grid->nz = nz;
  grid->subgrids = ptr_slist_new();
  GC_register_finalizer(grid, patch_grid_free, grid, NULL, NULL);
  return grid;
}

void patch_grid_add_subgrid(patch_grid_t* grid, patch_grid_t* subgrid)
{
  ASSERT(bbox_contains_bbox(&grid->domain, &subgrid->domain));
  ptr_slist_append(grid->subgrids, subgrid);
}

static bool subgrid_is_properly_nested(patch_grid_t* subgrid, 
                                       patch_grid_t* parent)
{
  real_t Lx = parent->domain.x2 - parent->domain.x1;
  real_t Ly = parent->domain.y2 - parent->domain.y1;
  real_t Lz = parent->domain.z2 - parent->domain.z1;
  real_t dx = Lx / parent->nx, dy = Ly / parent->ny, dz = Lz / parent->nz;
  bool nested = (fabs(dx * subgrid->domain.x1 - (int)(dx * subgrid->domain.x1)) < 1e-14);
  nested = nested && (fabs(dx * subgrid->domain.x2 - (int)(dx * subgrid->domain.x2)) < 1e-14);
  nested = nested && (fabs(dx * subgrid->domain.y1 - (int)(dy * subgrid->domain.y1)) < 1e-14);
  nested = nested && (fabs(dx * subgrid->domain.y2 - (int)(dy * subgrid->domain.y2)) < 1e-14);
  nested = nested && (fabs(dx * subgrid->domain.z1 - (int)(dz * subgrid->domain.z1)) < 1e-14);
  nested = nested && (fabs(dx * subgrid->domain.z2 - (int)(dz * subgrid->domain.z2)) < 1e-14);
  return nested;
}

bool patch_grid_is_properly_nested(patch_grid_t* grid)
{
  bool nested = true;
  ptr_slist_node_t* iter = grid->subgrids->front;
  if (iter == NULL)
    return true;
  else
  {
    while (iter != NULL)
    {
      nested = patch_grid_is_properly_nested(iter->value);
      if (nested)
        nested = subgrid_is_properly_nested(iter->value, grid);   
      if (!nested)
        return false;
      iter = iter->next;
    }
  }
  return nested;
}

