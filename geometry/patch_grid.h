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

#ifndef POLYMEC_PATCH_GRID_H
#define POLYMEC_GRID_H

#include "core/mesh.h"
#include "core/slist.h"

// This type describes a uniform Cartesian grid of a given resolution
// with a list of subgrids contained within it. Objects of this type are 
// garbage-collected.
typedef struct
{
  bbox_t domain;          // The domain spanned by this grid.
  int nx, ny, nz;         // Resolution in x, y, z.
  ptr_slist_t* subgrids;  // Grids within this one.
} patch_grid_t;

// Creates a new (garbage-collected) patch grid object.
patch_grid_t* patch_grid_new(bbox_t* domain, int nx, int ny, int nz);

// Adds a subgrid to the given patch grid.
void patch_grid_add_subgrid(patch_grid_t* grid, patch_grid_t* subgrid);

// Returns true if the patch grid contains subgrids that are properly nested, 
// false otherwise
bool patch_grid_is_properly_nested(patch_grid_t* grid);

#endif

