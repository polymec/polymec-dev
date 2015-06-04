// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/unordered_set.h"
#include "core/kd_tree.h"
#include "model/mesh_stencils.h"

static stencil_t* stencil_from_cells(const char* name, 
                                     int_array_t** stencil_cells,
                                     int num_cells,
                                     int num_ghost_cells,
                                     exchanger_t* ex)
{
  int* offsets = polymec_malloc(sizeof(int) * (num_cells+1));
  offsets[0] = 0;
  for (int i = 0; i < num_cells; ++i)
    offsets[i+1] = offsets[i] + stencil_cells[i]->size;
  int* indices = polymec_malloc(sizeof(int) * offsets[num_cells]);
  int k = 0;
  for (int i = 0; i < num_cells; ++i)
  {
    int_array_t* neighbors = stencil_cells[i];
    for (int i = 0; i < neighbors->size; ++i, ++k)
      indices[k] = neighbors->data[i];
    ASSERT(k == offsets[i+1]);
    int_array_free(neighbors);
  }
  polymec_free(stencil_cells);
  return unweighted_stencil_new(name, num_cells, offsets, indices, num_ghost_cells, ex);
}

stencil_t* cell_star_stencil_new(mesh_t* mesh, int radius)
{
  ASSERT(radius > 0);

  // First, we'll exchange the mesh cell centers to make sure they're 
  // consistent.
  exchanger_t* cell_ex = mesh_exchanger(mesh);
  exchanger_exchange(cell_ex, mesh->cell_centers, 3, 0, MPI_REAL_T);

  // First we will make a mapping from each cell to its list of 
  // neighboring cells.
  int_array_t** stencil_cells = polymec_malloc(sizeof(int_array_t*) * mesh->num_cells);
  for (int cell = 0; cell < mesh->num_cells; ++cell)
  {
    int_array_t* neighbors = int_array_new();
    int pos = 0, opp_cell;
    while (mesh_cell_next_neighbor(mesh, cell, &pos, &opp_cell))
    {
      if (opp_cell != -1)
        int_array_append(neighbors, opp_cell);
    }
    stencil_cells[cell] = neighbors;
  }

  // The exchanger for this stencil is the same as that for the mesh.
  exchanger_t* ex = exchanger_clone(mesh_exchanger(mesh));

  // Create the stencil from the sets.
  char name[1025];
  snprintf(name, 1024, "cell star stencil (R = %d)", radius);
  stencil_t* s = stencil_from_cells(name, stencil_cells, mesh->num_cells, mesh->num_ghost_cells, ex);

  // Now augment it.
  for (int r = 0; r < radius-1; ++r)
    stencil_augment(s);

  return s;
}

stencil_t* cell_halo_stencil_new(mesh_t* mesh, int radius)
{
  ASSERT(radius > 0);

  // Construct a slightly-bigger star stencil.
  stencil_t* s = cell_star_stencil_new(mesh, radius+1);

  // Rename it.
  polymec_free(s->name);
  char name[1025];
  snprintf(name, 1024, "cell halo stencil (R = %d)", radius);
  s->name = string_dup(name);

  // Now filter out the non-halo cells. 
  // FIXME

  return s;
}

