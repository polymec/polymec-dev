// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/unordered_set.h"
#include "core/kd_tree.h"
#include "model/polymesh_stencils.h"

static stencil_t* stencil_from_cells(const char* name,
                                     int_array_t** stencil_cells,
                                     int num_cells,
                                     int num_ghost_cells,
                                     exchanger_t* ex)
{
  int* offsets = polymec_malloc(sizeof(int) * (num_cells+1));
  offsets[0] = 0;
  for (int i = 0; i < num_cells; ++i)
    offsets[i+1] = offsets[i] + (int)stencil_cells[i]->size;
  int* indices = polymec_malloc(sizeof(int) * offsets[num_cells]);
  int k = 0;
  for (int i = 0; i < num_cells; ++i)
  {
    int_array_t* neighbors = stencil_cells[i];
    for (int j = 0; j < neighbors->size; ++j, ++k)
      indices[k] = neighbors->data[j];
    ASSERT(k == offsets[i+1]);
    int_array_free(neighbors);
  }
  polymec_free(stencil_cells);
  return stencil_new(name, num_cells, offsets, indices, num_ghost_cells, ex);
}

stencil_t* cell_star_stencil_new(polymesh_t* mesh, int radius)
{
  ASSERT(radius > 0);

  if (radius > 1)
    polymec_not_implemented("cell_star_stencil_new (radius > 1)");

  // First, we'll exchange the mesh cell centers to make sure they're
  // consistent.
  exchanger_t* cell_ex = polymesh_exchanger(mesh, POLYMESH_CELL);
  exchanger_exchange(cell_ex, mesh->cell_centers, 3, 0, MPI_REAL_T);

  // First we will make a mapping from each cell to its list of
  // neighboring cells.
  int_array_t** stencil_cells = polymec_malloc(sizeof(int_array_t*) * mesh->num_cells);
  for (int cell = 0; cell < mesh->num_cells; ++cell)
  {
    int_array_t* neighbors = int_array_new();
    int pos = 0, opp_cell;
    while (polymesh_cell_next_neighbor(mesh, cell, &pos, &opp_cell))
    {
      if (opp_cell != -1)
        int_array_append(neighbors, opp_cell);
    }
    stencil_cells[cell] = neighbors;
  }

  // The exchanger for this stencil is the same as that for the mesh.
  exchanger_t* ex = exchanger_clone(polymesh_exchanger(mesh, POLYMESH_CELL));

  // Create the stencil from the sets.
  char name[1025];
  snprintf(name, 1024, "cell star stencil (R = %d)", radius);
  stencil_t* s = stencil_from_cells(name, stencil_cells, mesh->num_cells, mesh->num_ghost_cells, ex);

//  // Now augment it.
//  for (int r = 0; r < radius-1; ++r)
//    stencil_augment(s);

  return s;
}

#if 0
noreturn stencil_t* cell_halo_stencil_new(polymesh_t* mesh)
{
  POLYMEC_NOT_IMPLEMENTED;

  // Construct a 2-deep star stencil to start from.
  stencil_t* star = cell_star_stencil_new(mesh, 1);

  // Now find the extra cells. The additional cells of a cell i are the
  // neighbors of the neighbors of i that are
  // (1) not connected to i
  // (2) connected to at least two of the neighbors of i.

  // First, make a list of the local neighbors we are adding for each cell.
  int_ptr_unordered_map_t* extra_neighbors = int_ptr_unordered_map_new();
  int_unordered_set_t* neighbors_of_i = int_unordered_set_new();
  for (int i = 0; i < star->num_indices; ++i)
  {
    // Make a list of the neighbors of i so we can filter them out.
    int posj = 0, j;
    while (stencil_next(star, i, &posj, &j, NULL))
      int_unordered_set_insert(neighbors_of_i, j);

    // Now traverse the neighbors of neighbors of i looking for extra halo
    // neighbors.
    posj = 0;
    while (stencil_next(star, i, &posj, &j, NULL))
    {
      // Populate the conn array with numbers of connections of each neighbor
      // of i's neighbors to neighbors of i (to evaluate (2) above).
      int num_n_o_n = stencil_size(star, j);
      int n_o_n[num_n_o_n];
      stencil_get_neighbors(star, j, n_o_n);
      int conn[num_n_o_n];
      memset(conn, 0, sizeof(int) * num_n_o_n);
      for (int kk = 0; kk < num_n_o_n; ++kk)
      {
        int k = n_o_n[kk];
        if (!int_unordered_set_contains(neighbors_of_i, k)) // (1)
        {
          int posl = 0, l;
          while (stencil_next(star, k, &posl, &l, NULL))
          {
            if (int_unordered_set_contains(neighbors_of_i, l))
              ++conn[kk];
          }
        }
      }

      // Now record all the neighbors of neighbors of i that are connected
      // to 2 or more neighbors of i. (2)
      for (int kk = 0; kk < num_n_o_n; ++kk)
      {
        if (conn[kk] >= 2)
        {
          int_array_t** halo_neighbors_p = (int_array_t**)int_ptr_unordered_map_get(extra_neighbors, i);
          int_array_t* halo_neighbors;
          if (halo_neighbors_p == NULL)
          {
            halo_neighbors = int_array_new();
            int_ptr_unordered_map_insert_with_v_dtor(extra_neighbors, i, halo_neighbors, DTOR(int_array_free));
          }
          else
            halo_neighbors = *halo_neighbors_p;
          int k = n_o_n[kk];
          int_array_append(halo_neighbors, k);
        }
      }
    }

    // Get ready to go again.
    int_unordered_set_clear(neighbors_of_i);
  }

  // Make sure that these extra neighbors are communicated between parallel
  // domains. First we traverse all send cells and make lists of extra send
  // cells to add. Then we communicate the numbers of additional sent cells
  // to the receiving processes. Finally we count up the additional number of
  // ghost cells and make sure they are reflected in our exchanger.

  // Traverse the send cells of the star stencil's exchanger and make lists
  // of extra send cells to add.

  // Communicate the numbers of additional sent cells to our neighboring
  // processes.
  int** send_data = (int**)exchanger_create_metadata_send_arrays(star->ex, MPI_INT, 1);
  int** receive_data = (int**)exchanger_create_metadata_send_arrays(star->ex, MPI_INT, 1);

  // Count up the needed ghost cells and make sure they are reflected in our
  // exchanger.
  int num_halo_ghosts = stencil_num_ghosts(star);
  // FIXME

  // Assemble the materials for the halo stencil.
  int* halo_offsets = NULL;
  int* halo_indices = NULL;
  exchanger_t* halo_ex = exchanger_new(exchanger_comm(star->ex));

  // Create the halo stencil.
  stencil_t* halo = stencil_new("cell halo stencil",
                                star->num_indices,
                                halo_offsets, halo_indices,
                                num_halo_ghosts, halo_ex);

  // Clean up.
  exchanger_free_metadata_send_arrays(star->ex, (void**)send_data);
  exchanger_free_metadata_receive_arrays(star->ex, (void**)receive_data);
  int_unordered_set_free(neighbors_of_i);
  int_ptr_unordered_map_free(extra_neighbors);
  stencil_free(star);

  return halo;
}
#endif

