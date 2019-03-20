// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/table.h"
#include "geometry/planar_polymesh.h"

// This function rounds the given number up to the nearest power of 2.
static int round_to_pow2(int x)
{
  int y = 2;
  while (y < x) y *= 2;
  return y;
}

planar_polymesh_t* planar_polymesh_new(int num_cells, int num_edges, int num_nodes)
{
  ASSERT(num_cells >= 0);
  ASSERT(num_edges >= 0);
  ASSERT(num_nodes >= 0);

  planar_polymesh_t* mesh = polymec_malloc(sizeof(planar_polymesh_t));

  // NOTE: We round stored elements up to the nearest power of 2.

  // Allocate cell information.
  mesh->num_cells = num_cells;
  mesh->cell_edge_offsets = polymec_calloc(num_cells+1, sizeof(int));
  mesh->cell_edge_cap = round_to_pow2(6 * num_cells);
  mesh->cell_edges = polymec_calloc(mesh->cell_edge_cap, sizeof(int));

  // Allocate edge information.
  mesh->num_edges = num_edges;
  mesh->edge_cells = polymec_malloc(2*num_edges*sizeof(int));
  for (int e = 0; e < 2*mesh->num_edges; ++e)
    mesh->edge_cells[e] = -1;
  mesh->edge_nodes = polymec_calloc(2*num_edges, sizeof(int));

  // Allocate node information.
  mesh->num_nodes = num_nodes;
  mesh->nodes = polymec_calloc(num_nodes, sizeof(point2_t));

  // Allocate tagging mechanisms.
  mesh->cell_tags = tagger_new();
  mesh->edge_tags = tagger_new();
  mesh->node_tags = tagger_new();

  return mesh;
}

planar_polymesh_t* planar_polymesh_new_with_cell_type(int num_cells,
                                                      int num_edges,
                                                      int num_nodes,
                                                      int num_edges_per_cell)
{
  planar_polymesh_t* mesh = planar_polymesh_new(num_cells, num_edges, num_nodes);

  // Set up connectivity metadata.
  for (int c = 0; c < mesh->num_cells+1; ++c)
    mesh->cell_edge_offsets[c] = num_edges_per_cell*c;
  planar_polymesh_reserve_connectivity_storage(mesh);
  return mesh;
}

void planar_polymesh_free(planar_polymesh_t* mesh)
{
  ASSERT(mesh != NULL);

  // Destroy tags.
  tagger_free(mesh->node_tags);
  tagger_free(mesh->edge_tags);
  tagger_free(mesh->cell_tags);

  // Destroy nodes.
  polymec_free(mesh->nodes);

  // Destroy connectivity.
  polymec_free(mesh->edge_nodes);
  polymec_free(mesh->edge_cells);
  polymec_free(mesh->cell_edges);
  polymec_free(mesh->cell_edge_offsets);

  // Destroy the mesh itself.
  polymec_free(mesh);
}

bool planar_polymesh_is_valid(planar_polymesh_t* mesh, char** reason)
{
  static char _reason[1025];

  // All cells must have at least 3 edges.
  for (int c = 0; c < mesh->num_cells; ++c)
  {
    if (planar_polymesh_cell_num_edges(mesh, c) < 3)
    {
      snprintf(_reason, 1024, "Polygonal cell %d has only %d edges.",
               c, planar_polymesh_cell_num_edges(mesh, c));
      if (reason != NULL)
        *reason = _reason;
      return false;
    }
  }

  // Check cell-edge topology.
  for (int c = 0; c < mesh->num_cells; ++c)
  {
    int pos = 0, e;
    while (planar_polymesh_cell_next_edge(mesh, c, &pos, &e))
    {
      if (e >= mesh->num_edges)
      {
        snprintf(_reason, 1024, "Cell %d has invalid edge %d (mesh has only "
                 "%d edges).", c, e, mesh->num_edges);
        if (reason != NULL)
          *reason = _reason;
        return false;
      }
      if ((mesh->edge_cells[2*e] != c) && (mesh->edge_cells[2*e+1] != c))
      {
        snprintf(_reason, 1024, "Cell %d has edge %d in its list of edges, "
                 "but that edge does not have that cell in its list.", c, e);
        if (reason != NULL)
          *reason = _reason;
        return false;
      }
    }
  }
  for (int e = 0; e < mesh->num_edges; ++e)
  {
    int cell1 = mesh->edge_cells[2*e];
    if (cell1 == -1)
    {
      snprintf(_reason, 1024, "Edge %d has no cells in its list.", e);
      if (reason != NULL)
        *reason = _reason;
      return false;
    }
    if (cell1 >= mesh->num_cells)
    {
      snprintf(_reason, 1024, "Edge %d has invalid cell %d (mesh has only %d "
               "cells).", e, cell1, mesh->num_cells);
      if (reason != NULL)
        *reason = _reason;
      return false;
    }

    int pos = 0, ee;
    bool found_edge = false;
    while (planar_polymesh_cell_next_edge(mesh, cell1, &pos, &ee))
    {
      if (ee == e)
      {
        found_edge = true;
        break;
      }
    }
    if (!found_edge)
    {
      snprintf(_reason, 1024, "Edge %d has cell %d in its list of cells, but "
              "that cell does not have that edge in its list.", e, cell1);
      if (reason != NULL)
        *reason = _reason;
      return false;
    }

    int cell2 = mesh->edge_cells[2*e+1];
    if (cell2 != -1)
    {
      if (cell2 >= mesh->num_cells)
      {
        snprintf(_reason, 1024, "Edge %d has invalid cell %d (mesh has only "
                 "%d cells).", e, cell2, mesh->num_cells);
        if (reason != NULL)
          *reason = _reason;
        return false;
      }

      pos = 0;
      found_edge = false;
      while (planar_polymesh_cell_next_edge(mesh, cell2, &pos, &ee))
      {
        if (ee == e)
        {
          found_edge = true;
          break;
        }
      }
      if (!found_edge)
      {
        snprintf(_reason, 1024, "Edge %d has cell %d in its list of cells, but "
                "that cell does not have that edge in its list.", e, cell2);
        if (reason != NULL)
          *reason = _reason;
        return false;
      }
    }
  }

  // Check edge-node topology.
  for (int e = 0; e < mesh->num_edges; ++e)
  {
    int n1 = mesh->edge_nodes[2*e];
    int n2 = mesh->edge_nodes[2*e+1];
    if (n1 == n2)
    {
      snprintf(_reason, 1024, "Edge %d has the same node (%d) at both ends!",
               e, n1);
      if (reason != NULL)
        *reason = _reason;
      return false;
    }
    else if ((n1 < 0) || (n1 >= mesh->num_nodes))
    {
      snprintf(_reason, 1024, "Edge %d has invalid first node %d (mesh has "
               "only %d nodes).", e, n1, mesh->num_nodes);
      if (reason != NULL)
        *reason = _reason;
      return false;
    }
    else if ((n2 < 0) || (n2 >= mesh->num_nodes))
    {
      snprintf(_reason, 1024, "Edge %d has invalid second node %d (mesh has "
               "only %d nodes).", e, n2, mesh->num_nodes);
      if (reason != NULL)
        *reason = _reason;
      return false;
    }
  }

  return true;
}

planar_polymesh_t* planar_polymesh_clone(planar_polymesh_t* mesh)
{
  planar_polymesh_t* clone = planar_polymesh_new(mesh->num_cells,
                                                 mesh->num_edges,
                                                 mesh->num_nodes);

  // Connectivity metadata.
  memcpy(clone->cell_edge_offsets, mesh->cell_edge_offsets, sizeof(int)*(mesh->num_cells+1));
  planar_polymesh_reserve_connectivity_storage(clone);

  // Actual connectivity.
  int num_cell_edges = clone->cell_edge_offsets[clone->num_cells];
  memcpy(clone->cell_edges, mesh->cell_edges, sizeof(int)*num_cell_edges);
  memcpy(clone->edge_cells, mesh->edge_cells, sizeof(int)*2*clone->num_edges);
  memcpy(clone->edge_nodes, mesh->edge_nodes, sizeof(int)*2*clone->num_edges);

  // Node stuff.
  memcpy(clone->nodes, mesh->nodes, sizeof(point2_t)*clone->num_nodes);

  return clone;
}

void planar_polymesh_reserve_connectivity_storage(planar_polymesh_t* mesh)
{
  // Make sure metadata is in order.
  int num_cell_edges = mesh->cell_edge_offsets[mesh->num_cells];
  ASSERT(num_cell_edges >= 3*mesh->num_cells);

  if (mesh->cell_edge_cap <= num_cell_edges)
  {
    while (mesh->cell_edge_cap <= num_cell_edges)
      mesh->cell_edge_cap *= 2;
    mesh->cell_edges = polymec_realloc(mesh->cell_edges, sizeof(int) * mesh->cell_edge_cap);
  }
}

void planar_polymesh_cell_get_nodes(planar_polymesh_t* mesh,
                                    int cell,
                                    int* nodes)
{
  ASSERT(cell >= 0);
  ASSERT(cell < mesh->num_cells);
  int n = 0;
  for (int e = mesh->cell_edge_offsets[cell]; e < mesh->cell_edge_offsets[cell+1]; ++e, ++n)
  {
    int edge = mesh->cell_edges[e];
    if (edge >= 0)
      nodes[n] = mesh->edge_nodes[2*edge];
    else
      nodes[n] = mesh->edge_nodes[2*(~edge)+1];
  }
}

polygon_t* planar_polymesh_cell_polygon(planar_polymesh_t* mesh, int cell)
{
  ASSERT(cell >= 0);
  ASSERT(cell < mesh->num_cells);
  size_t num_vertices = mesh->cell_edge_offsets[cell+1] -
                        mesh->cell_edge_offsets[cell];
  int nodes[num_vertices];
  planar_polymesh_cell_get_nodes(mesh, cell, nodes);
  point2_t vertices[num_vertices];
  for (size_t v = 0; v < num_vertices; ++v)
    vertices[v] = mesh->nodes[nodes[v]];
  return polygon_new(vertices, num_vertices);
}

void planar_polymesh_cell_get_polygon(planar_polymesh_t* mesh,
                                      int cell,
                                      polygon_t* polygon)
{
  ASSERT(cell >= 0);
  ASSERT(cell < mesh->num_cells);
  size_t num_vertices = mesh->cell_edge_offsets[cell+1] -
                        mesh->cell_edge_offsets[cell];
  int nodes[num_vertices];
  planar_polymesh_cell_get_nodes(mesh, cell, nodes);
  point2_t vertices[num_vertices];
  for (size_t v = 0; v < num_vertices; ++v)
    vertices[v] = mesh->nodes[nodes[v]];
  polygon_set_vertices(polygon, vertices, num_vertices);
}

adj_graph_t* graph_from_planar_polymesh_cells(planar_polymesh_t* mesh)
{
  // Create a graph whose vertices are the mesh's cells.
  adj_graph_t* g = adj_graph_new(MPI_COMM_SELF, mesh->num_cells);

  // Allocate space in the graph for the edges (edges connecting cells).
  for (int i = 0; i < mesh->num_cells; ++i)
  {
    // How many edges don't have opposite cells?
    int outer_edges = 0;
    for (int j = mesh->cell_edge_offsets[i]; j < mesh->cell_edge_offsets[i+1]; ++j)
    {
      int e = mesh->cell_edges[j];
      if (e < 0)
        e = ~e;
      if (mesh->edge_cells[2*e+1] == -1)
        ++outer_edges;
    }
    adj_graph_set_num_edges(g, i, mesh->cell_edge_offsets[i+1] - mesh->cell_edge_offsets[i] - outer_edges);
  }

  // Now fill in the edges.
  for (int i = 0; i < mesh->num_cells; ++i)
  {
    int* edges = adj_graph_edges(g, i);
    int offset = 0;
    for (int j = mesh->cell_edge_offsets[i]; j < mesh->cell_edge_offsets[i+1]; ++j)
    {
      int e = mesh->cell_edges[j];
      if (e < 0)
        e = ~e;
      if (mesh->edge_cells[2*e+1] != -1)
      {
        int c = (i == mesh->edge_cells[2*e]) ? mesh->edge_cells[2*e+1] : mesh->edge_cells[2*e];
        edges[offset] = c;
        ++offset;
      }
    }
  }

  return g;
}

