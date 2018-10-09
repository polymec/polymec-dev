// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

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
  hex_map_t* hexes = hexagon(radius);

  // Traverse the cells and pick off elements.
  int_array_t* cell_edges = int_array_new_with_size(6*hexes->size);
  for (size_t i = 0; i < cell_edges->size; ++i)
    cell_edges->data[i] = -1;
  int_array_t* edge_cells = int_array_new();
  int_array_t* edge_nodes = int_array_new();
  point2_array_t* nodes = point2_array_new();
  int pos = 0, cell_index;
  hex_t* hex;
  while (hex_map_next(hexes, &pos, &hex, &cell_index))
  {
    // Add edges and nodes for this cell.
    for (int dir = 0; dir < 6; ++dir)
    {
      if (cell_edges->data[6*cell_index+dir] == -1) // no edge/node yet
      {
        // Fetch the neighboring hex in this direction.
        hex_t n;
        hex_map_get_neighbor(hexes, hex, dir, &n);
        int neighbor_index = *hex_map_get(hexes, &n);
        int neighbor_dir = (dir + 3) % 6;

        // Add the edge to both hexes.
        int edge_index = (int)(edge_cells->size/2);
        cell_edges->data[6*cell_index+dir] = edge_index;
        cell_edges->data[6*neighbor_index+neighbor_dir] = edge_index;

        // Add both hexes to the edge.
        int_array_append(edge_cells, cell_index);
        int_array_append(edge_cells, neighbor_index);

        // Create the edge's first node.
        if (dir == 0)
        {
          int_array_append(edge_nodes, nodes->size);
          point2_t n1;
          hex_get_node_position(hex, dir, &n1);
          point2_array_append(nodes, n1);
        }
        else
        {
          // This edge's first node is the second node of the edge in 
          // the prior direction.
          int prior_dir = (dir + 5) % 6;
          int prior_edge = cell_edges->data[6*cell_index+prior_dir];
          int n1 = edge_nodes->data[2*prior_edge+1];
          int_array_append(edge_nodes, n1);
        }

        // Create the edge's second node (if needed).
        if (dir < 5)
        {
          int_array_append(edge_nodes, nodes->size);
          int next_dir = (dir+1)%6;
          point2_t n2;
          hex_get_node_position(hex, next_dir, &n2);
          point2_array_append(nodes, n2);
        }
        else
        {
          // The last node in the last edge of the cell is the first 
          // node in the cell's first edge.
          int first_edge = cell_edges->data[6*cell_index];
          int n1 = edge_nodes->data[2*first_edge];
          int_array_append(edge_nodes, n1);
        }
      }
    }
  }

  // Create our planar polymesh.
  planar_polymesh_t* mesh = planar_polymesh_new_with_cell_type((int)hexes->size,
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
  hex_map_free(hexes);

  return mesh;
}

