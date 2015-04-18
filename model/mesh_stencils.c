// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/unordered_set.h"
#include "core/hilbert.h"
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

stencil_t* cell_star_stencil_new(mesh_t* mesh, int radius)
{
  ASSERT(radius > 0);

  // First, we'll exchange the mesh cell centers to make sure they're 
  // consistent.
  exchanger_t* cell_ex = mesh_exchanger(mesh);
  exchanger_exchange(cell_ex, mesh->cell_centers, 3, 0, MPI_REAL_T);

  // Now we stick all of our cell centers into a kd-tree.
  kd_tree_t* cell_tree = kd_tree_new(mesh->cell_centers, mesh->num_cells);

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
  {
    POLYMEC_NOT_IMPLEMENTED;
#if 0
    // We count up our non-ghost neighbors and exchange with our neighbors.
    exchanger_t* cell_ex = mesh_exchanger(mesh);
    int num_interior_neighbors[mesh->num_cells];
    int max_interior_neighbors = 0;
    bbox_t bbox = {.x1 = FLT_MAX, .x2 = -FLT_MAX,
                   .y1 = FLT_MAX, .y2 = -FLT_MAX,
                   .z1 = FLT_MAX, .z2 = -FLT_MAX};
    for (int cell = 0; cell < mesh->num_cells; ++cell)
    {
      num_interior_neighbors[cell] = 0;
      bbox_grow(&bbox, &mesh->cell_centers[cell]);
      for (int i = 0; i < neighbors[cell]->size; ++i)
      {
        int neighbor = neighbors[cell]->data[2*i];
        if (neighbor < mesh->num_cells)
        {
          ++num_interior_neighbors[cell];
          if (num_interior_neighbors[cell] > max_interior_neighbors)
            max_interior_neighbors = num_interior_neighbors[cell];
        }
      }
    }
    exchanger_exchange(num_interior_neighbors, 1, 0, MPI_INT);
    MPI_Allreduce(MPI_IN_PLACE, &max_interior_neighbors, 1, MPI_INT, MPI_MAX, mesh->comm);

    // We'll need a Hilbert curve for sorting cell indices objectively.
    hilbert_t* hilbert = hilbert_new(&bbox);

    // Now pack the actual data for exchanging. 
    int* interior_neighbors = polymec_malloc(sizeof(int) * 2 * max_interior_neighbors * mesh->num_cells);
    for (int cell = 0; cell < mesh->num_cells; ++cell)
    {
      int i = 0;
      point_t centers[num_interior_neighbors[cell]];
      for (int i = 0; i < neighbors[cell]->size; ++i)
      {
        int neighbor = neighbors[cell]->data[2*i];
      int pos = 0, pos1 = 0, neighbor, distance;
      while (int_unordered_set_next(neighbors, &pos, &neighbor))
      {
        if (neighbor < mesh->num_cells)
        {
          interior_neighbors[2*(max_interior_neighbors*cell+i)] = neighbor;
          centers[i] = mesh->cell_centers[interior_neighbors[max_interior_neighbors*cell+i]];
          ++i;
        }
      }

      // Sort the interior neighbors into Hilbert order.
      hilbert_sort(hilbert, centers, 
                   &interior_neighbors[max_interior_neighbors*cell],
                   num_interior_neighbors[cell]);
    }

    // Now exchange.
    exchanger_exchange(interior_neighbors, max_interior_neighbors, 0, MPI_INT);

    // Now construct the exchanger.
    int_ptr_unordered_map_t* send_map = int_ptr_unordered_map_new();
    int_ptr_unordered_map_t* receive_map = int_ptr_unordered_map_new();
    for (int cell = 0; cell < mesh->num_cells; ++cell)
    {
      // Does this cell have any ghost neighbors?
      if (num_interior_neighbors[cell] < stencil_sets[cell]->size)
      {
        // Pick through the neighbors for this cell and 
      }
    }
    exchanger_set_sends(ex, send_map);
    exchanger_set_sends(ex, receive_map);

    // Clean up.
    int_ptr_unordered_map_free(send_map);
    int_ptr_unordered_map_free(receive_map);
    polymec_free(interior_neighbors);
#endif
  }

  // Clean up.
  kd_tree_free(cell_tree);

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

  // Now we stick all of our cell centers into a kd-tree.
  kd_tree_t* cell_tree = kd_tree_new(mesh->cell_centers, mesh->num_cells);

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
  exchanger_t* ex;
  if (radius == 1)
  {
    // The exchanger for this stencil is the same as that for the mesh.
    ex = exchanger_clone(mesh_exchanger(mesh));
  }
  else
  {
    POLYMEC_NOT_IMPLEMENTED;
  }

  // Create the stencil from the sets.
  char name[1025];
  snprintf(name, 1024, "cell halo stencil (R = %d)", radius);
  return stencil_from_cells(name, stencil_cells, mesh->num_cells, ex);
}

