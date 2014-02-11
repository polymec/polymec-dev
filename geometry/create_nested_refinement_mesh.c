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
#include "geometry/create_nested_refinement_mesh.h"

static void nested_grid_free(void* ctx, void* dummy)
{
  nested_grid_t* grid = ctx;
  ptr_slist_free(grid->subgrids);
}

nested_grid_t* nested_grid_new(bbox_t* domain, int nx, int ny, int nz)
{
  ASSERT(domain->x1 < domain->x2);
  ASSERT(domain->y1 < domain->y2);
  ASSERT(domain->z1 < domain->z2);
  ASSERT(nx > 0);
  ASSERT(ny > 0);
  ASSERT(nz > 0);
  nested_grid_t* grid = GC_MALLOC(sizeof(nested_grid_t));
  grid->domain = *domain;
  grid->nx = nx;
  grid->ny = ny;
  grid->nz = nz;
  grid->subgrids = ptr_slist_new();
  GC_register_finalizer(grid, nested_grid_free, grid, NULL, NULL);
  return grid;
}

void nested_grid_add_subgrid(nested_grid_t* grid, nested_grid_t* subgrid)
{
  ASSERT(bbox_contains_bbox(&grid->domain, &subgrid->domain));
  ptr_slist_append(grid->subgrids, subgrid);
}

mesh_t* create_nested_refinement_mesh(MPI_Comm comm, 
                                      nested_grid_t* coarse_grid)
{
  return NULL;
}
 
