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

#include "model/stencil.h"
#include "core/unordered_set.h"
#include "core/array.h"

stencil_t* stencil_new(const char* name, int num_indices, 
                       int* offsets, int* indices, real_t* weights,
                       exchanger_t* ex)
{
  ASSERT(num_indices > 0);
  ASSERT(offsets != NULL);
  ASSERT(indices != NULL);
  stencil_t* s = polymec_malloc(sizeof(stencil_t));
  s->name = string_dup(name);
  s->num_indices = num_indices;
  s->offsets = offsets;
  s->indices = indices;
  s->weights = weights;
  s->ex = ex;
  return s;
}

stencil_t* unweighted_stencil_new(const char* name, int num_indices, 
                                  int* offsets, int* indices, exchanger_t* ex)
{
  return stencil_new(name, num_indices, offsets, indices, NULL, ex);
}

void stencil_free(stencil_t* stencil)
{
  polymec_free(stencil->name);
  polymec_free(stencil->offsets);
  polymec_free(stencil->indices);
  if (stencil->weights != NULL)
    polymec_free(stencil->weights);
  exchanger_free(stencil->ex);
  polymec_free(stencil);
}

stencil_t* cell_face_stencil_new(mesh_t* mesh)
{
  // Read off the stencil from the mesh.
  // FIXME: This is broken at the moment, since we have -1's in there.
  int* offsets = polymec_malloc(sizeof(int) * (mesh->num_cells+1));
  memcpy(offsets, mesh->cell_face_offsets, sizeof(int) * (mesh->num_cells+1));
  int* indices = polymec_malloc(sizeof(int) * offsets[mesh->num_cells]);
  for (int i = 0; i < offsets[mesh->num_cells]; ++i)
    indices[i] = mesh_face_opp_cell(mesh, mesh->cell_faces[offsets[i]], i);
  exchanger_t* ex = exchanger_clone(mesh_exchanger(mesh));
  return unweighted_stencil_new("cell-face stencil", mesh->num_cells, offsets, 
                                indices, ex);
}

stencil_t* cell_edge_stencil_new(mesh_t* mesh)
{
  // First we will make a big clumsy mapping from each cell to its list of 
  // neighboring cells.
  int_array_t** stencil_map = polymec_malloc(sizeof(int_array_t*) * mesh->num_cells);
  int_unordered_set_t* cell_edges = int_unordered_set_new();
  int_unordered_set_t* cell_neighbors = int_unordered_set_new();
  exchanger_t* ex = exchanger_new(mesh->comm);
  for (int cell = 0; cell < mesh->num_cells; ++cell)
  {
    int_array_t* neighbors = int_array_new();

    // Make a list of all the edges in the cell.
    int_unordered_set_clear(cell_edges);
    int pos = 0, face;
    while (mesh_cell_next_face(mesh, cell, &pos, &face))
    {
      int fpos = 0, edge;
      while (mesh_face_next_edge(mesh, face, &fpos, &edge))
        int_unordered_set_insert(cell_edges, edge);
    }
//printf("cell %d has edges ", cell);
//int epos = 0, e;
//while (int_unordered_set_next(cell_edges, &epos, &e))
//printf("%d ", e);
//printf("\n");

    int opp_cell;
    pos = 0;
    int_unordered_set_clear(cell_neighbors);
    while (mesh_cell_next_neighbor(mesh, cell, &pos, &opp_cell))
    {
      if (opp_cell != -1)
      {
        // This one's definitely a neighbor.
        if (!int_unordered_set_contains(cell_neighbors, opp_cell))
        {
          int_array_append(neighbors, opp_cell);
          int_unordered_set_insert(cell_neighbors, opp_cell);
        }
        
        // Each neighbor of this opposite cell that shares an edge with our 
        // original cell is a neighbor of that cell.
        int opos = 0, face;
        while (mesh_cell_next_face(mesh, opp_cell, &opos, &face))
        {
          int ncell = mesh_face_opp_cell(mesh, face, opp_cell);
          if ((ncell != -1) && (ncell != cell) && 
              !int_unordered_set_contains(cell_neighbors, ncell))
          {
            int fpos = 0, edge;
            while (mesh_face_next_edge(mesh, face, &fpos, &edge))
            {
//printf("%d: %d\n", ncell, edge);
              if (int_unordered_set_contains(cell_edges, edge) &&
                  !int_unordered_set_contains(cell_neighbors, ncell))
              {
//printf("%d: has edge in common with %d (via face %d)\n", cell, ncell, face);
                int_array_append(neighbors, ncell);
                int_unordered_set_insert(cell_neighbors, ncell);
              }
            }
          }
        }
      }
    }

    stencil_map[cell] = neighbors;
  }
  int_unordered_set_free(cell_edges);
  int_unordered_set_free(cell_neighbors);

  // FIXME: Still have to construct the exchanger for this stencil.

  // Create the stencil from the map.
  int* offsets = polymec_malloc(sizeof(int) * (mesh->num_cells+1));
  offsets[0] = 0;
  for (int i = 0; i < mesh->num_cells; ++i)
    offsets[i+1] = offsets[i] + stencil_map[i]->size;
  int* indices = polymec_malloc(sizeof(int) * offsets[mesh->num_cells]);
  int k = 0;
  for (int i = 0; i < mesh->num_cells; ++i)
  {
    int_array_t* neighbors = stencil_map[i];
    for (int j = 0; j < neighbors->size; ++j, ++k)
      indices[k] = neighbors->data[j];
    ASSERT(k == offsets[i+1]);
    int_array_free(neighbors);
  }
  polymec_free(stencil_map);
  return unweighted_stencil_new("cell-edge stencil", mesh->num_cells, offsets, 
                                indices, ex);
}

stencil_t* cell_node_stencil_new(mesh_t* mesh)
{
  // First we will make a big clumsy mapping from each cell to its list of 
  // neighboring cells.
  int_array_t** stencil_map = polymec_malloc(sizeof(int_array_t*) * mesh->num_cells);
  int_unordered_set_t* cell_nodes = int_unordered_set_new();
  int_unordered_set_t* cell_neighbors = int_unordered_set_new();
  int_unordered_set_t* cell_neighbor_candidates = int_unordered_set_new();
  exchanger_t* ex = exchanger_new(mesh->comm);
  for (int cell = 0; cell < mesh->num_cells; ++cell)
  {
    int_array_t* neighbors = int_array_new();

    // Make a list of all the nodes in the cell.
    int_unordered_set_clear(cell_nodes);
    int pos = 0, face;
    while (mesh_cell_next_face(mesh, cell, &pos, &face))
    {
      int fpos = 0, node;
      while (mesh_face_next_node(mesh, face, &fpos, &node))
        int_unordered_set_insert(cell_nodes, node);
    }
//printf("cell %d has nodes ", cell);
//int npos = 0, n;
//while (int_unordered_set_next(cell_nodes, &npos, &n))
//printf("%d ", n);
//printf("\n");

    // Now gather neighbor candidates.
    int opp_cell;
    pos = 0;
    int_unordered_set_clear(cell_neighbor_candidates);
    while (mesh_cell_next_neighbor(mesh, cell, &pos, &opp_cell))
    {
      if (opp_cell != -1)
      {
        int_unordered_set_insert(cell_neighbor_candidates, opp_cell);
        
        // Each neighbor of this opposite cell is also a candidate.
        int opos = 0, ncell;
        while (mesh_cell_next_neighbor(mesh, opp_cell, &opos, &ncell))
        {
          if ((ncell != -1) && (ncell != cell))
          {
            int_unordered_set_insert(cell_neighbor_candidates, ncell);
          
            // Finally, each neighbor of this opposite cell is also a candidate.
            int npos = 0, nncell;
            while (mesh_cell_next_neighbor(mesh, ncell, &npos, &nncell))
            {
              if ((nncell != -1) && (nncell != cell))
                int_unordered_set_insert(cell_neighbor_candidates, nncell);
            }
          }
        }
      }
    }

    // Now search through the candidates and determine which ones share nodes
    // with the original cell.
    int cpos = 0, ncell;
    int_unordered_set_clear(cell_neighbors);
    while (int_unordered_set_next(cell_neighbor_candidates, &cpos, &ncell))
    {
      int fpos = 0, face;
      bool cell_is_neighbor = false;
      while (mesh_cell_next_face(mesh, ncell, &fpos, &face))
      {
        int npos = 0, node;
        while (mesh_face_next_node(mesh, face, &npos, &node))
        {
          if (int_unordered_set_contains(cell_nodes, node) && 
              !int_unordered_set_contains(cell_neighbors, ncell))
          {
            int_array_append(neighbors, ncell);
            int_unordered_set_insert(cell_neighbors, ncell);
            cell_is_neighbor = true;
            break;
          }
        }
        if (cell_is_neighbor) break;
      }
    }

    stencil_map[cell] = neighbors;
  }
  int_unordered_set_free(cell_nodes);
  int_unordered_set_free(cell_neighbors);
  int_unordered_set_free(cell_neighbor_candidates);

  // FIXME: Still have to construct the exchanger for this stencil.

  // Create the stencil from the map.
  int* offsets = polymec_malloc(sizeof(int) * (mesh->num_cells+1));
  offsets[0] = 0;
  for (int i = 0; i < mesh->num_cells; ++i)
    offsets[i+1] = offsets[i] + stencil_map[i]->size;
  int* indices = polymec_malloc(sizeof(int) * offsets[mesh->num_cells]);
  int k = 0;
  for (int i = 0; i < mesh->num_cells; ++i)
  {
    int_array_t* neighbors = stencil_map[i];
    for (int j = 0; j < neighbors->size; ++j, ++k)
      indices[k] = neighbors->data[j];
    ASSERT(k == offsets[i+1]);
    int_array_free(neighbors);
  }
  polymec_free(stencil_map);
  return unweighted_stencil_new("cell-node stencil", mesh->num_cells, offsets, 
                                indices, ex);
}

