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
                                     exchanger_t* ex)
{
  int* offsets = polymec_malloc(sizeof(int) * (num_cells+1));
  offsets[0] = 0;
  for (int i = 0; i < num_cells; ++i)
    offsets[i+1] = offsets[i] + stencil_cells[i]->size/2;
  int* indices = polymec_malloc(sizeof(int) * offsets[num_cells]);
  int k = 0;
  for (int i = 0; i < num_cells; ++i)
  {
    int_array_t* neighbors = stencil_cells[i];
    ASSERT((neighbors->size % 2) == 0);
    for (int i = 0; i < neighbors->size/2; ++i, ++k)
      indices[k] = neighbors->data[2*i];
    ASSERT(k == offsets[i+1]);
    int_array_free(neighbors);
  }
  polymec_free(stencil_cells);
  return unweighted_stencil_new(name, num_cells, offsets, indices, ex);
}

static exchanger_t* exchanger_from_cells(mesh_t* mesh,
                                         int_array_t** stencil_cells, 
                                         int radius)
{
  int num_cells = mesh->num_cells;
  exchanger_t* cell_ex = mesh_exchanger(mesh);

  // Exchange the number of interior neighbors with other processes.
  int total_num_cells = num_cells + mesh->num_ghost_cells;
  int our_num_interior_neighbors[total_num_cells];
  int max_interior_neighbors = 0;
  for (int i = 0; i < num_cells; ++i)
  {
    our_num_interior_neighbors[i] = 0;
    for (int j = 0; j < stencil_cells[j]->size; ++j)
    {
      int neighbor = stencil_cells[j]->data[2*i];
      if (neighbor < mesh->num_cells)
        ++our_num_interior_neighbors[i];
      if (our_num_interior_neighbors[i] > max_interior_neighbors)
        max_interior_neighbors = our_num_interior_neighbors[i];
    }
  }
  exchanger_exchange(cell_ex, our_num_interior_neighbors, 1, 0, MPI_INT);
  MPI_Allreduce(MPI_IN_PLACE, &max_interior_neighbors, 1, MPI_INT, MPI_MAX, mesh->comm);

  // Now pack the actual data for exchanging. 
  int* remote_distances = polymec_malloc(sizeof(int) * max_interior_neighbors * total_num_cells);
  memset(remote_distances, 0, sizeof(int) * max_interior_neighbors * total_num_cells);
  for (int i = 0; i < num_cells; ++i)
  {
    int k = 0;
    for (int j = 0; j < stencil_cells[i]->size; ++j)
    {
      int neighbor = stencil_cells[i]->data[2*j];
      int distance = stencil_cells[i]->data[2*j+1];
      if (neighbor >= mesh->num_cells)
      {
        remote_distances[max_interior_neighbors*i+k] = distance;
        ++k;
      }
    }
  }

  // Now exchange.
  exchanger_exchange(cell_ex, remote_distances, max_interior_neighbors, 0, MPI_INT);

  // Now construct the send/receive maps for the exchanger.
  int_ptr_unordered_map_t* send_map = int_ptr_unordered_map_new();
  int_ptr_unordered_map_t* receive_map = int_ptr_unordered_map_new();
  int next_ghost_index = mesh->num_ghost_cells;
  for (int i = 0; i < mesh->num_cells; ++i)
  {
    // Pick through the neighbors for this cell and add the 
    // ones that are close enough.
    for (int j = 0; j < stencil_cells[i]->size; ++j)
    {
      int neighbor = stencil_cells[i]->data[2*j];
      int base_distance = stencil_cells[i]->data[2*j+1];
      if (neighbor < num_cells)
      {
        // We set up our send cells here.
      }
      else
      {
        int remote_distance = remote_distances[max_interior_neighbors*i+j];
        if ((remote_distance > 0) && // remote neighbor!
            (base_distance + remote_distance <= radius))
        {
          // We set up our receive cells here.
        }
      }
    }
  }

  // Create the exchanger.
  exchanger_t* ex = exchanger_new(mesh->comm);
  exchanger_set_sends(ex, send_map);
  exchanger_set_receives(ex, receive_map);

  // Clean up.
  int_ptr_unordered_map_free(send_map);
  int_ptr_unordered_map_free(receive_map);
  polymec_free(remote_distances);

  return ex;
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
  int_unordered_set_t* neighbor_set = int_unordered_set_new();
  for (int cell = 0; cell < mesh->num_cells; ++cell)
  {
    int_array_t* neighbors = int_array_new();

    // Start at this cell and journey outward through faces.
    int distance = 0;
    int_array_t* cells = int_array_new();
    int_array_t* new_cells = int_array_new();
    int_array_append(cells, cell);
    while (distance < radius)
    {
      for (int i = 0; i < cells->size; ++i)
      {
        int pos = 0, opp_cell;
        while (mesh_cell_next_neighbor(mesh, cells->data[i], &pos, &opp_cell))
        {
          if ((opp_cell != -1) && 
              (!int_unordered_set_contains(neighbor_set, opp_cell)) && 
              (opp_cell != cell))
          {
            int_unordered_set_insert(neighbor_set, opp_cell);
            int_array_append(neighbors, opp_cell);
            int_array_append(neighbors, distance + 1);
            int_array_append(new_cells, opp_cell);
          }
        }
      }

      ++distance;
      int_array_free(cells);
      cells = new_cells;
    }
    stencil_cells[cell] = neighbors;
    int_unordered_set_clear(neighbor_set);
  }
  int_unordered_set_free(neighbor_set);

  // Construct the exchanger for this stencil.
  exchanger_t* ex;
  if (radius == 1)
  {
    // The exchanger for this stencil is the same as that for the mesh.
    ex = exchanger_clone(mesh_exchanger(mesh));
  }
  else
    ex = exchanger_from_cells(mesh, stencil_cells, radius);

  // Create the stencil from the sets.
  char name[1025];
  snprintf(name, 1024, "cell star stencil (R = %d)", radius);
  return stencil_from_cells(name, stencil_cells, mesh->num_cells, ex);
}

stencil_t* cell_halo_stencil_new(mesh_t* mesh, int radius)
{
  ASSERT(radius > 0);

  // First, we'll exchange the mesh cell centers to make sure they're 
  // consistent.
  exchanger_t* cell_ex = mesh_exchanger(mesh);
  exchanger_exchange(cell_ex, mesh->cell_centers, 3, 0, MPI_REAL_T);

  // First we will make a mapping from each cell to its list of 
  // neighboring cells.
  int_array_t** stencil_cells = polymec_malloc(sizeof(int_array_t*) * mesh->num_cells);
  int_unordered_set_t* neighbor_set = int_unordered_set_new();
  int_unordered_set_t* node_set = int_unordered_set_new();
  for (int cell = 0; cell < mesh->num_cells; ++cell)
  {
    int_array_t* neighbors = int_array_new();

    // Make a list of all the nodes attached to this cell.
    int fpos = 0, face;
    while (mesh_cell_next_face(mesh, cell, &fpos, &face))
    {
      int npos = 0, node;
      while (mesh_face_next_node(mesh, face, &npos, &node))
        int_unordered_set_insert(node_set, node);
    }

    // Start at this cell and journey outward through faces.
    int distance = 0;
    int_array_t* cells = int_array_new();
    int_array_t* new_cells = int_array_new();
    int_array_append(cells, cell);
    while (distance <= radius)
    {
      for (int i = 0; i < cells->size; ++i)
      {
        int pos = 0, opp_cell;
        while (mesh_cell_next_neighbor(mesh, cells->data[i], &pos, &opp_cell))
        {
          if ((opp_cell != -1) && 
              (!int_unordered_set_contains(neighbor_set, opp_cell)) &&
              (opp_cell != cell))
          {
            int_unordered_set_insert(neighbor_set, opp_cell);

            // Does this cell share a node with our original cell?
            bool shares_node = false;
            int fpos = 0, face;
            while (mesh_cell_next_face(mesh, opp_cell, &fpos, &face))
            {
              int npos = 0, node;
              while (mesh_face_next_node(mesh, face, &npos, &node))
              {
                if (int_unordered_set_contains(node_set, node))
                  shares_node = true;
              }
            }

            if (shares_node)
            {
              int_array_append(neighbors, opp_cell);
              int_array_append(neighbors, distance);
              int_array_append(new_cells, opp_cell);
            }
            else if ((distance + 1) <= radius)
            {
              int_array_append(neighbors, opp_cell);
              int_array_append(neighbors, distance + 1);
              int_array_append(new_cells, opp_cell);
            }
          }
        }
      }

      ++distance;
      int_array_free(cells);
      cells = new_cells;
    }
    stencil_cells[cell] = neighbors;
    int_unordered_set_clear(neighbor_set);
    int_unordered_set_clear(node_set);
  }
  int_unordered_set_free(neighbor_set);
  int_unordered_set_free(node_set);

  // Construct the exchanger for this stencil.
  exchanger_t* ex = exchanger_from_cells(mesh, stencil_cells, radius);

  // Create the stencil from the sets.
  char name[1025];
  snprintf(name, 1024, "cell halo stencil (R = %d)", radius);
  return stencil_from_cells(name, stencil_cells, mesh->num_cells, ex);
}

