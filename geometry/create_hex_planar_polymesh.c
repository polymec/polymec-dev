// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/array.h"
#include "core/unordered_map.h"
#include "geometry/create_hex_planar_polymesh.h"

// Hex coordinate type (using axial coordinates as defined at
// https://www.redblobgames.com/grids/hexagons).
typedef struct
{
  int q, r;
} hex_t;

static hex_t* hex_new(int q, int r)
{
  hex_t* hex = polymec_malloc(sizeof(hex_t));
  hex->q = q;
  hex->r = r;
  return hex;
}

static void hex_free(hex_t* hex)
{
  polymec_free(hex);
}

static inline int hex_hash(hex_t* hex)
{
  return int_pair_hash((int*)hex);
}

static inline bool hex_equals(hex_t* x, hex_t* y)
{
  return int_pair_equals((int*)x, (int*)y);
}

DEFINE_UNORDERED_MAP(hex_map, hex_t*, int, hex_hash, hex_equals)

// Create a hexagonal array of, um, hexagons, with the given radius.
static hex_map_t* hexagon(int radius)
{
  hex_map_t* hexes = hex_map_new();
  int i = 0;
  for (int q = -radius; q <= radius; ++q)
  {
    int r1 = MAX(-radius, -q - radius);
    int r2 = MIN(radius, -q + radius);
    for (int r = r1; r <= r2; ++r, ++i)
      hex_map_insert_with_k_dtor(hexes, hex_new(q, r), i, hex_free);
  }
  return hexes;
}

// Starting angles for x- and y-aligned hexes.
// Constructs the neighboring hex for the given hex, in the given direction.
static inline void hex_map_get_neighbor(hex_map_t* hexes, 
                                        hex_t* hex,
                                        int direction, 
                                        hex_t* neighbor)
{
  // q, r displacements by direction
  static int dq[6] = {1, 1, 0, -1, -1, 0};
  static int dr[6] = {0, -1, -1, 0, 1, 1};
  neighbor->q = hex->q + dq[direction]; 
  neighbor->r = hex->r + dr[direction];
}

static void hex_get_node_position(hex_t* hex, 
                                  int direction, 
                                  real_t h,
                                  point2_t* x)
{
  // Find the hex's center.
  const real_t sqrt3 = sqrt(3.0);
  real_t x0 = h * 1.5 * hex->q;
  real_t y0 = -h * (0.5*sqrt3*hex->q + sqrt3*hex->r);

  // Get to the node from there.
  real_t theta = (direction-1) * M_PI/3.0;
  x->x = x0 + h*cos(theta);
  x->y = y0 + h*sin(theta);
printf("NODE POSITION FOR (%d, %d), dir %d: (%g, %g)\n", hex->q, hex->r, direction, x->x, x->y);
}

planar_polymesh_t* create_hex_planar_polymesh(size_t radius, 
                                              real_t h)
{
  ASSERT(h > 0.0);

  // Build a hexagonal array of hexes.
  hex_map_t* hex_map = hexagon((int)radius);

  // Traverse the cells and pick off elements.
  hex_t hex_inv_map[hex_map->size];
  int_array_t* cell_edges = int_array_new_with_size(6*hex_map->size);
  int no_edge = INT_MAX;
  for (size_t i = 0; i < cell_edges->size; ++i)
    cell_edges->data[i] = no_edge;
  int_array_t* edge_cells = int_array_new();
  int pos = 0, cell_index;
  hex_t* hex;
  while (hex_map_next(hex_map, &pos, &hex, &cell_index))
  {
    // Add edges for this cell.
    for (int dir = 0; dir < 6; ++dir)
    {
      if (cell_edges->data[6*cell_index+dir] == no_edge) // no edge yet
      {
        // Fetch the neighboring hex in this direction.
        hex_t n;
        hex_map_get_neighbor(hex_map, hex, dir, &n);
        int* neighbor_index_p = hex_map_get(hex_map, &n);
        int neighbor_index = -1;
        if (neighbor_index_p != NULL) // not on a boundary
          neighbor_index = *neighbor_index_p;

        // Add the (oriented) edge to both hexes.
        int edge_index = (int)(edge_cells->size/2);
        cell_edges->data[6*cell_index+dir] = edge_index;
        if (neighbor_index != -1)
        {
          int neighbor_dir = (dir + 3) % 6;
          cell_edges->data[6*neighbor_index+neighbor_dir] = ~edge_index;
        }

        // Add both hexes to the edge.
        int_array_append(edge_cells, cell_index);
        int_array_append(edge_cells, neighbor_index);
      }
    }

    // Construct the inverse mapping.
    hex_inv_map[cell_index] = *hex;
  }

  // Now create nodes for each of the edges.
  size_t num_edges = edge_cells->size / 2;
  int_array_t* edge_nodes = int_array_new_with_size(2*num_edges);
  for (size_t e = 0; e < 2*num_edges; ++e)
    edge_nodes->data[e] = -1;
  point2_array_t* nodes = point2_array_new();
  pos = 0;
  while (hex_map_next(hex_map, &pos, &hex, &cell_index))
  {
printf("Cell %d (%d, %d)\n", cell_index, hex->q, hex->r);
    for (int dir = 0; dir < 6; ++dir)
    {
      int edge_index = cell_edges->data[6*cell_index+dir];
      int edge_node1_index, edge_node2_index;
      if (edge_index < 0)
      {
        edge_index = ~edge_index;
        edge_node1_index = 2*edge_index+1;
        edge_node2_index = 2*edge_index;
      }
      else
      {
        edge_node1_index = 2*edge_index;
        edge_node2_index = 2*edge_index+1;
      }
printf(" edge #%d\n", dir);

      // Access the neighbor cell and determine the edge indices of the 
      // nodes incident on either side. 
      int neighbor_cell = edge_cells->data[2*edge_index+1];
      int neighbor_edge1_index = no_edge, neighbor_edge2_index = no_edge;
      int neighbor_edge_node1_index = -1, 
          neighbor_edge_node2_index = -1;
      if (neighbor_cell != -1)
      {
        hex_t* nhex = &(hex_inv_map[neighbor_cell]);
printf(" Found neighbor cell (%d, %d)\n", nhex->q, nhex->r);
        int opp_dir = (dir + 3) % 6;

        // Find the edge incident to node 1 within the neighbor cell.
        neighbor_edge1_index = cell_edges->data[6*neighbor_cell+(opp_dir+1)%6];
        if (neighbor_edge1_index < 0)
        {
          neighbor_edge1_index = ~neighbor_edge1_index;
          neighbor_edge_node1_index = 2*neighbor_edge1_index+1;
        }
        else
          neighbor_edge_node1_index = 2*neighbor_edge1_index;

        // Find the edge incident to node 2 within the neighbor cell.
        neighbor_edge2_index = cell_edges->data[6*neighbor_cell+(opp_dir+5)%6];
        if (neighbor_edge2_index < 0)
        {
          neighbor_edge2_index = ~neighbor_edge2_index;
          neighbor_edge_node2_index = 2*neighbor_edge2_index;
        }
        else
          neighbor_edge_node2_index = 2*neighbor_edge2_index+1;
      }

      // Find the first node for this edge.
      if (edge_nodes->data[edge_node1_index] == -1) // no first node yet
      {
        // Create the first node for this edge.
        int n1 = (int)nodes->size;
        point2_t x1;
        hex_get_node_position(hex, dir, h, &x1);
        point2_array_append(nodes, x1);
printf("Created node %d (#1) for edge %d\n", n1, edge_index);

        // Hook up the node to this edge.
        edge_nodes->data[edge_node1_index] = n1;

        // This node is the second node of the previous edge in this cell.
        int prev_edge_index = cell_edges->data[6*cell_index+(dir+5)%6];
        int prev_edge_node2_index;
        if (prev_edge_index < 0)
        {
          prev_edge_index = ~prev_edge_index;
          prev_edge_node2_index = 2*prev_edge_index;
        }
        else
          prev_edge_node2_index = 2*prev_edge_index+1;
printf("Attaching node %d as #2 for previous edge %d\n", n1, prev_edge_index);
        ASSERT(edge_nodes->data[prev_edge_node2_index] == -1);
        edge_nodes->data[prev_edge_node2_index] = n1;
      }
      else printf("Edge %d already had node #1: %d\n", edge_index, edge_nodes->data[edge_node1_index]);

      // If this edge connects this cell to a neighbor cell, hook up the node 
      // to the other incident edge in that neighbor cell.
      if ((neighbor_cell != -1) && (edge_nodes->data[neighbor_edge1_index] == -1))
{
printf("Attaching node %d to neighbor edge %d\n", edge_nodes->data[edge_node1_index], neighbor_edge1_index);
        edge_nodes->data[neighbor_edge_node1_index] = edge_nodes->data[edge_node1_index];
}

      // Find the second node for this edge (if needed).
      if (edge_nodes->data[edge_node2_index] == -1) // no second node yet
      {
        // Create the second node for this edge.
        int n2 = (int)nodes->size;
        point2_t x2;
        hex_get_node_position(hex, (dir+1)%6, h, &x2);
        point2_array_append(nodes, x2);
printf("Created node %d (#2) for edge %d\n", n2, edge_index);

        // Hook up the node to this edge.
        edge_nodes->data[edge_node2_index] = n2;

        // This node is the first node in the next edge of this cell.
        int next_edge_index = cell_edges->data[6*cell_index+(dir+1)%6];
        int next_edge_node1_index;
        if (next_edge_index < 0)
        {
          next_edge_index = ~next_edge_index;
          next_edge_node1_index = 2*next_edge_index+1;
        }
        else
          next_edge_node1_index = 2*next_edge_index;
printf("Attaching node %d as #1 for next edge %d\n", n2, next_edge_index);
        ASSERT(edge_nodes->data[next_edge_node1_index] == -1);
        edge_nodes->data[next_edge_node1_index] = n2;
      }
else printf("Edge %d already had node #2: %d\n", edge_index, edge_nodes->data[edge_node2_index]);

      // Attach the neighbor's incident edge to this cell.
      if ((neighbor_cell != -1) && (edge_nodes->data[neighbor_edge_node2_index] == -1))
{
printf("Attaching node %d to neighbor edge %d\n", edge_nodes->data[edge_node2_index], neighbor_edge2_index);
        edge_nodes->data[neighbor_edge_node2_index] = edge_nodes->data[edge_node2_index];
}
    }
  }

  // Create our planar polymesh.
  planar_polymesh_t* mesh = planar_polymesh_new_with_cell_type((int)hex_map->size,
                                                               (int)num_edges,
                                                               (int)nodes->size, 
                                                               6);
  memcpy(mesh->cell_edges, cell_edges->data, sizeof(int) * cell_edges->size);
  memcpy(mesh->edge_nodes, edge_nodes->data, sizeof(int) * edge_nodes->size);
  memcpy(mesh->edge_cells, edge_cells->data, sizeof(int) * edge_cells->size);
  memcpy(mesh->nodes, nodes->data, sizeof(point2_t) * nodes->size);

  // Clean up.
  int_array_free(cell_edges);
  int_array_free(edge_nodes);
  int_array_free(edge_cells);
  point2_array_free(nodes);
  hex_map_free(hex_map);

  return mesh;
}

