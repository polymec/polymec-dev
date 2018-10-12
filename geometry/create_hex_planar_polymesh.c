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
DEFINE_ARRAY(hex_array, hex_t*)

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

static void hex_get_node_position(hex_t* hex, int direction, point2_t* x)
{
//static const real_t start_angles[2] = {0.0,       // x-aligned
//                                       M_PI/2.0}; // y-aligned
  // FIXME
}

planar_polymesh_t* create_hex_planar_polymesh(hex_lattice_align_t alignment,
                                              size_t radius, real_t h)
{
  ASSERT(h > 0.0);

  // Build a hexagonal array of hexes.
  hex_map_t* hex_map = hexagon((int)radius);

  // Traverse the cells and pick off elements.
  hex_t hex_inv_map[hex_map->size];
  int_array_t* cell_edges = int_array_new_with_size(6*hex_map->size);
  for (size_t i = 0; i < cell_edges->size; ++i)
    cell_edges->data[i] = -1;
  int_array_t* edge_cells = int_array_new();
  int pos = 0, cell_index;
  hex_t* hex;
  while (hex_map_next(hex_map, &pos, &hex, &cell_index))
  {
    // Add edges for this cell.
    for (int dir = 0; dir < 6; ++dir)
    {
      if (cell_edges->data[6*cell_index+dir] == -1) // no edge/node yet
      {
        // Fetch the neighboring hex in this direction.
        hex_t n;
        hex_map_get_neighbor(hex_map, hex, dir, &n);
        int* neighbor_index_p = hex_map_get(hex_map, &n);
        int neighbor_index = -1;
        if (neighbor_index_p != NULL) // not on a boundary
          neighbor_index = *neighbor_index_p;

        int neighbor_dir = (dir + 3) % 6;

        // Add the (oriented) edge to both hexes.
        int edge_index = (int)(edge_cells->size/2);
        cell_edges->data[6*cell_index+dir] = edge_index;
        if (neighbor_index != -1)
          cell_edges->data[6*neighbor_index+neighbor_dir] = ~edge_index;

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
    for (int dir = 0; dir < 6; ++dir)
    {
      int edge_index = cell_edges->data[6*cell_index+dir];
      if (edge_index >= 0)
      {
        if (edge_nodes->data[2*edge_index] == -1)
        {
          // Create the first node for this edge.
          int n1 = (int)nodes->size;
          point2_t x1;
          hex_get_node_position(hex, dir, &x1);
          point2_array_append(nodes, x1);

          // Hook up the node to this edge.
          edge_nodes->data[2*edge_index] = n1;

          // Hook up the node to the other edge it belongs to in this cell.
          int prev_edge_index = cell_edges->data[6*cell_index+(dir+5)%6];
          ASSERT(edge_nodes->data[2*prev_edge_index+1] == -1);
          edge_nodes->data[2*prev_edge_index+1] = n1;
        }

        // Create the second node for this edge (if needed).
        if (edge_nodes->data[2*edge_index+1] == -1)
        {
          int n2 = (int)nodes->size;
          point2_t x2;
          hex_get_node_position(hex, (dir+1)%6, &x2);
          point2_array_append(nodes, x2);

          // Hook up the node to this edge.
          edge_nodes->data[2*edge_index+1] = n2;

          // Hook up the node to the other edge it belongs to in this cell.
          int next_edge_index = cell_edges->data[6*cell_index+(dir+1)%6];
          ASSERT(edge_nodes->data[2*next_edge_index+1] == -1);
          edge_nodes->data[2*next_edge_index+1] = n2;
        }
      }
    }
  }

  // Create our planar polymesh.
  planar_polymesh_t* mesh = planar_polymesh_new_with_cell_type((int)hex_map->size,
                                                               (int)cell_edges->size,
                                                               (int)nodes->size, 
                                                               6);
  memcpy(mesh->cell_edges, cell_edges->data, cell_edges->size);
  memcpy(mesh->edge_nodes, edge_nodes->data, edge_nodes->size);
  memcpy(mesh->edge_cells, edge_cells->data, edge_cells->size);
  memcpy(mesh->nodes, nodes->data, sizeof(point2_t) * nodes->size);

  // Clean up.
  int_array_free(cell_edges);
  int_array_free(edge_nodes);
  int_array_free(edge_cells);
  point2_array_free(nodes);
  hex_map_free(hex_map);

  return mesh;
}

