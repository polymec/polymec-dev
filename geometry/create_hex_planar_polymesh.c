// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "geometry/create_hex_planar_polymesh.h"

// Hex coordinate type (using axial coordinates as defined at
// https://www.redblobgames.com/grids/hexagons).
typedef struct
{
  int q, r;
} hex_t;

// Hex arithmetic.

// c1 + c2 -> sum
static void hex_add(const hex_t* c1, const hex_t* c2, hex_t* sum)
{
  sum->q = c1->q + c2->q;
  sum->r = c1->r + c2->r;
}

// c1 - c2 -> diff
//static void hex_sub(const hex_t* c1, const hex_t* c2, hex_t* diff)
//{
//  diff->q = c1->q - c2->q;
//  diff->r = c1->r - c2->r;
//}

// k * c -> prod
static void hex_scale(const hex_t* hex, int k, hex_t* prod)
{
  prod->q = k * hex->q;
  prod->r = k * hex->r;
}

// Distance between two hex cells c1 and c2.
static inline int hex_distance(hex_t* c1, hex_t* c2)
{
  return (ABS(c1->q - c2->q) + 
          ABS(c1->q + c1->r - c2->q - c2->r) + 
          ABS(c1->r - c2->r)) / 2;
}

// Direction vectors.
static const hex_t hex_directions[6] = 
  {{+1, 0}, {+1, -1}, {0, -1}, {-1, 0}, {-1, +1}, {0, +1}};

//static inline void hex_get_dir(int dir, hex_t* dir_vector)
//{
//  return hex_directions[dir];
//}

static inline void hex_get_neighbor(hex_t* hex, 
                                    int dir,
                                    hex_t* neighbor)
{
  hex_add(hex, &(hex_directions[dir]), neighbor);
}

#if 0
// Node spacings for alignments/directions, to be scaled by h.
static const point2_t node_spacings[2][6] = {
  // x-aligned hex lattice
  {{}},
  // y-aligned hex lattice
  {{}}
};
#endif

planar_polymesh_t* create_hex_planar_polymesh(hex_lattice_align_t alignment,
                                              size_t radius, real_t h)
{
  ASSERT(h > 0.0);

  // Traverse the cells, starting at the center and spiraling outward. 
  // Add edges and nodes as we go.
  int ncells = 0;
  int_array_t* cell_edges = int_array_new();
  int_array_t* edge_cells = int_array_new();
  int_array_t* edge_nodes = int_array_new();
  point2_array_t* nodes = point2_array_new();
  if (radius == 0)
  {
    // The mesh contains a single hex.
    ncells = 1;
    point2_t n[2][6]; // FIXME
    for (int e = 0; e < 6; ++e)
    {
      int_array_append(cell_edges, e);
      int_array_append(edge_nodes, e);
      int_array_append(edge_nodes, (e+1)%6);
      int_array_append(edge_cells, 0);
      int_array_append(edge_cells, -1);
      point2_array_append(nodes, n[(int)alignment][e]);
    }
  }
  else
  {
    // This traversal follows 
    // https://www.redblobgames.com/grids/hexagons/#rings.
    hex_t center = {0, 0};
    for (int rad = 0; rad <= radius; ++rad)
    {
      // move outward along direction 4
      hex_t hex, dr;
      hex_scale(&hex_directions[4], rad, &dr);
      hex_add(&center, &dr, &hex);

      // Traverse the ring of hexes at this radius.
      for (int dir = 0; dir < 6; ++dir)
      {
        for (int phi = 0; phi <= radius; ++phi, ++ncells)
        {
          // FIXME: Do stuff!
          // Move to the next hex in the ring.
          hex_get_neighbor(&hex, dir, &hex);
        }
      }
    }
  }

  // Create our planar polymesh.
  planar_polymesh_t* mesh = planar_polymesh_new_with_cell_type(ncells,
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

  return mesh;
}

