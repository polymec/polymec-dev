// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

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

void stencil_exchange(stencil_t* stencil, void* data, int stride, int tag, MPI_Datatype type)
{
  exchanger_exchange(stencil->ex, data, stride, tag, type);
}

int stencil_start_exchange(stencil_t* stencil, void* data, int stride, int tag, MPI_Datatype type)
{
  return exchanger_start_exchange(stencil->ex, data, stride, tag, type);
}

void stencil_finish_exchange(stencil_t* stencil, int token)
{
  exchanger_finish_exchange(stencil->ex, token);
}

exchanger_t* stencil_exchanger(stencil_t* stencil)
{
  return stencil->ex;
}

static size_t stencil_byte_size(void* obj)
{
  stencil_t* stencil = obj;
  
  // Data.
  size_t basic_storage = sizeof(int) + sizeof(char) * strlen(stencil->name) + 
                         sizeof(int) * (1 + stencil->num_indices + 1 + stencil->offsets[stencil->num_indices] + 1);
  if (stencil->weights != NULL)
    basic_storage += sizeof(real_t) * stencil->offsets[stencil->num_indices];
  
  // Exchanger-related storage.
  serializer_t* ex_s = exchanger_serializer();
  size_t ex_storage = serializer_size(ex_s, stencil->ex);
  ex_s = NULL;

  return basic_storage + ex_storage;
}

static void* stencil_byte_read(byte_array_t* bytes, size_t* offset)
{
  // Read the name.
  int name_len; 
  byte_array_read_ints(bytes, 1, &name_len, offset);
  char name[name_len+1];
  byte_array_read_chars(bytes, name_len, name, offset);
  name[name_len] = '\0';

  // Read the offsets, indices, weights.
  int num_indices, max_offset, num_weights;
  byte_array_read_ints(bytes, 1, &num_indices, offset);
  int* offsets = polymec_malloc(sizeof(int) * (num_indices+1));
  byte_array_read_ints(bytes, num_indices+1, offsets, offset);
  byte_array_read_ints(bytes, 1, &max_offset, offset);
  int* indices = polymec_malloc(sizeof(int) * max_offset);
  byte_array_read_ints(bytes, max_offset, indices, offset);
  byte_array_read_ints(bytes, 1, &num_weights, offset);
  real_t* weights = NULL;
  if (num_weights > 0)
    byte_array_read_real_ts(bytes, num_weights, weights, offset);

  // Exchanger stuff.
  serializer_t* ser = exchanger_serializer();
  exchanger_t* ex = serializer_read(ser, bytes, offset);

  return stencil_new(name, num_indices, offsets, indices, weights, ex);
}

static void stencil_byte_write(void* obj, byte_array_t* bytes, size_t* offset)
{
  stencil_t* stencil = obj;

  // Write the name.
  int name_len = strlen(stencil->name);
  byte_array_write_ints(bytes, 1, &name_len, offset);
  byte_array_write_chars(bytes, name_len, stencil->name, offset);

  // Write the offsets, indices, weights.
  byte_array_write_ints(bytes, 1, &stencil->num_indices, offset);
  byte_array_write_ints(bytes, stencil->num_indices+1, stencil->offsets, offset);
  int max_offset = stencil->offsets[stencil->num_indices];
  byte_array_write_ints(bytes, 1, &max_offset, offset);
  byte_array_write_ints(bytes, max_offset, stencil->indices, offset);
  if (stencil->weights != NULL)
  {
    byte_array_write_ints(bytes, 1, &max_offset, offset);
    byte_array_write_real_ts(bytes, max_offset, stencil->weights, offset);
  }
  else
  {
    int zero = 0;
    byte_array_write_ints(bytes, 1, &zero, offset);
  }

  // Exchanger.
  serializer_t* ser = exchanger_serializer();
  serializer_write(ser, stencil->ex, bytes, offset);
  ser = NULL;
}

serializer_t* stencil_serializer()
{
  return serializer_new("stencil", stencil_byte_size, stencil_byte_read, stencil_byte_write, DTOR(stencil_free));
}

stencil_t* cell_face_stencil_new(mesh_t* mesh)
{
  // First we will make a big clumsy mapping from each cell to its list of 
  // neighboring cells.
  int_array_t** stencil_map = polymec_malloc(sizeof(int_array_t*) * mesh->num_cells);
  for (int cell = 0; cell < mesh->num_cells; ++cell)
  {
    int_array_t* neighbors = int_array_new();

    int pos = 0, opp_cell;
    while (mesh_cell_next_neighbor(mesh, cell, &pos, &opp_cell))
    {
      if (opp_cell != -1)
        int_array_append(neighbors, opp_cell);
    }
    stencil_map[cell] = neighbors;
  }

  // The exchanger for this stencil is the same as that for the mesh.
  exchanger_t* ex = exchanger_clone(mesh_exchanger(mesh));

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
  return unweighted_stencil_new("cell-face stencil", mesh->num_cells, offsets, 
                                indices, ex);
}

// Here's a struct holding cell connectivity properties.
typedef struct
{
  int num_faces, num_edges, num_nodes;
  int faces[64], edges[128], nodes[256];
} cell_properties_t;

static stencil_t* create_layered_mesh_stencil(mesh_t* mesh,
                                              int num_layers,
                                              bool (*cell_matches)(void* context, 
                                                                   cell_properties_t* cell, 
                                                                   int cell_index),
                                              void* context)
{
  ASSERT(num_layers > 0);

  // The first layer is identical to the cell-face stencil.
  int layer = 1;
  stencil_t* s = cell_face_stencil_new(mesh);

  // Make a list of cells 

  // Send the 

  // Now we determine the cells that belong within the rest of the layers of 
  // the stencil.
  while (layer < num_layers)
  {
    ++layer;
  }

  // Send back all the cells that are attached to the cells near ghosts.

  return s;
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

adj_graph_t* graph_from_point_cloud_and_stencil(point_cloud_t* points, 
                                                stencil_t* stencil)
{
  ASSERT(stencil_num_indices(stencil) == points->num_points);

  // Create a graph whose vertices are the cloud's points.
  int rank, nproc;
  MPI_Comm_size(points->comm, &nproc);
  MPI_Comm_rank(points->comm, &rank);
  adj_graph_t* g = adj_graph_new(points->comm, points->num_points);

  // Allocate space in the graph for the edges (stencil size).
  int num_points = points->num_points;
  for (int i = 0; i < num_points; ++i)
    adj_graph_set_num_edges(g, i, stencil_size(stencil, i));

  // Now fill in the edges.
  for (int i = 0; i < num_points; ++i)
  {
    int* edges = adj_graph_edges(g, i);
    int offset = 0, pos = 0, j;
    while (stencil_next(stencil, i, &pos, &j, NULL))
      edges[offset++] = j;
  }

  return g;
}

