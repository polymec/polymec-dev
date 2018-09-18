// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "geometry/create_hex_planar_polymesh.h"

typedef struct 
{
  // The alignment of the hex lattice.
  hex_lattice_align_t alignment;
  // The radius of the hex lattice (in outward cells from the origin).
  size_t radius;

  // Bookkeeping stuff. Look away!
  size_t nc;
  int* nc_for_r;
} hex_lattice_t;

static hex_lattice_t* hex_lattice_new(hex_lattice_align_t alignment,
                                      size_t radius)
{
  hex_lattice_t* l = polymec_malloc(sizeof(hex_lattice_t));
  l->alignment = alignment;
  l->radius = radius;

  // Cache numbers of cells stored within our given radius.
  l->nc = 1;
  l->nc_for_r = polymec_malloc(sizeof(int) * (radius+1));
  for (size_t r = 0; r < radius; ++r)
  {
    l->nc += 6*(r+1);
    l->nc_for_r[r] = (int)(l->nc);
  }

  return l;
}

static void hex_lattice_free(hex_lattice_t* l)
{
  polymec_free(l->nc_for_r);
  polymec_free(l);
}

static inline size_t hex_lattice_num_cells(hex_lattice_t* l)
{
  return l->nc;
}

static inline size_t hex_lattice_num_edges(hex_lattice_t* l)
{
  size_t n = 6;
  for (size_t r = 0; r < l->radius; ++r)
    n += 6*(3+2*r); 
  return n;
}

static inline size_t hex_lattice_num_nodes(hex_lattice_t* l)
{
  // Same as the number of edges!
  return hex_lattice_num_edges(l);
}

static inline void hex_lattice_get_displacement(hex_lattice_t* l, 
                                                int dir, int* dq, int* dr)
{
  ASSERT(dir >= 0);
  ASSERT(dir < 6);

  // Taken from https://www.redblobgames.com/grids/hexagons/#neighbors-axial
  static const int q_disps[6] = {+1, +1,  0, -1, -1,  0};
  static const int r_disps[6] = { 0, -1, -1,  0, +1, +1};
  *dq = q_disps[dir];
  *dr = r_disps[dir];
}

static inline void hex_lattice_cell_get_neighbor(hex_lattice_t* l, 
                                                 int q, int r, int dir,
                                                 int* q1, int* r1) 
{
  hex_lattice_get_displacement(l, dir, q1, r1);
  *q1 += q;
  *r1 += r;
}

static inline int hex_lattice_cell_distance(hex_lattice_t* l, 
                                            int q1, int r1,
                                            int q2, int r2)
{
  return (ABS(q1 - q2) + ABS(q1 + r1 - q2 - r2) + ABS(r1 - r2)) / 2;
}

static inline int hex_lattice_cell_wedge(hex_lattice_t* l, int q, int r)
{
  // FIXME: Need to figure out the relationship between q, r, s, and x, y, z.
  int s = -q - r;
  int q_turf = q - r;
  int r_turf = r - s;
  int s_turf = s - q;
  int max = MAX(ABS(q_turf), MAX(ABS(r_turf), ABS(s_turf)));
  return (max == q_turf) ? (q_turf > 0) ? 0 : 3
                         : (max == r_turf) ? (r_turf > 0) ? 1 : 4
                                           : (s_turf > 0) ? 2 : 5;
}

static int hex_lattice_cell(hex_lattice_t* l, int q, int r)
{
  // Find the distance of this point from the origin.
  int dist = hex_lattice_cell_distance(l, 0, 0, q, r);
  if (dist == 0)
    return 0;

  // Unit radial displacement "vector".
  static int rad_hat[2];
  hex_lattice_get_displacement(l, 0, &rad_hat[0], &rad_hat[1]);

  // Increment (q1, r1) till we reach (q, r), starting at the first position 
  // at the right distance.
  int q1 = dist * rad_hat[0];
  int r1 = dist * rad_hat[1];
  int index = l->nc_for_r[dist];
  if ((q1 != q) || (r1 != r))
  {
    for (int d = 0; d < 6; ++d)
    {
      int dir = (d + 2) % 6;
      for (int i = 0; i < dist; ++i)
      {
        hex_lattice_cell_get_neighbor(l, q1, r1, dir, &q1, &r1);
        ++index;
        if ((q1 == q) && (r1 == r))
          break;
      }
    }
  }

  return index; 
}

#if 0
static void hex_lattice_get_cell_pair(hex_lattice_t* l, 
                                      int index, int* q, int* r)
{
  *q = 0; 
  *r = 0;
  if (index == 0)
    return;

  // Find the radius of the cell.
  int radius = 0;
  while (l->nc_for_r[radius+1] < index)
    ++radius;

  // Update q and r to reflect the radius.
  static int rad_hat[2];
  hex_lattice_get_displacement(l, 0, &rad_hat[0], &rad_hat[1]);
  *q = radius * rad_hat[0];
  *r = radius * rad_hat[1];

  // Now increment (q, r) till we reach the desired index.
  // at the right distance.
  int current_index = l->nc_for_r[radius];
  if (current_index != index)
  {
    for (int d = 0; d < 6; ++d)
    {
      int dir = (d + 2) % 6;
      for (int i = 0; i < radius; ++i)
      {
        hex_lattice_cell_get_neighbor(l, *q, *r, dir, q, r);
        ++current_index;
        if (current_index == index)
          break;
      }
    }
  }
}
#endif

static int hex_lattice_cell_edge(hex_lattice_t* l, int q, int r, int dir)
{
  // Find the distance of this point from the origin.
  int dist = hex_lattice_cell_distance(l, 0, 0, q, r);
  if (dist == 0)
    return dir;

  // Find the wedge that the cell is in.
  int wedge = hex_lattice_cell_wedge(l, q, r);

  // Find the edge by traversing our way there.
  // FIXME
  return dist*wedge + dir;
}

static inline int hex_lattice_cell_node(hex_lattice_t* l, int q, int r, int dir) 
{
  return hex_lattice_cell_edge(l, q, r, dir);
}

static bool hex_lattice_next_cell(hex_lattice_t* l, int* pos, int* q, int* r)
{
  if (*pos >= l->nc)
    return false;
  else 
  {
    if (*pos == 0)
      *q = *r = 0;
    else
    {
      // (Inspired by https://www.redblobgames.com/grids/hexagons/#rings)

      // Figure out what ring we're in.
      int cell = *pos;
      int radius = 0;
      while (cell >= l->nc_for_r[radius])
        ++radius;

      // Are we graduating to the next ring?
      if ((cell + 1) >= l->nc_for_r[radius])
      {
        // Move one cell "outward."
        int dir = 0; 
        hex_lattice_cell_get_neighbor(l, *q, *r, dir, q, r);
      }
      else
      {
        // We're staying in this ring, so figure out the direction we need 
        // to go next.
        int ring_offset = cell + 1 - l->nc_for_r[radius];
        int dir = ring_offset / radius + 1;

        // Update q and r.
        hex_lattice_cell_get_neighbor(l, *q, *r, dir, q, r);
      }
    }
    ++(*pos);
    return true;
  }
}

planar_polymesh_t* create_hex_planar_polymesh(hex_lattice_align_t alignment,
                                              size_t radius, real_t h)
{
  ASSERT(h > 0.0);

  // Create a hex lattice that we use to read off the connectivity.
  hex_lattice_t* lattice = hex_lattice_new(alignment, radius);
  size_t num_cells = hex_lattice_num_cells(lattice);
  size_t num_edges = hex_lattice_num_edges(lattice);
  size_t num_nodes = hex_lattice_num_nodes(lattice);

  // Create our planar polymesh.
  planar_polymesh_t* mesh = planar_polymesh_new_with_cell_type((int)num_cells,
                                                               (int)num_edges,
                                                               (int)num_nodes, 
                                                               6);

  // Traverse the cells, starting at (0, 0) and spiraling outward. 
  // Read off edges and nodes as we go.
  int cpos = 0, q, r;
  while (hex_lattice_next_cell(lattice, &cpos, &q, &r))
  {
    int cell = hex_lattice_cell(lattice, q, r);

    // Go around the cell along the +/- q/r/s directions and collect 
    // its edges and nodes.
    for (int dir = 0; dir < 6; ++dir)
    {
      // Edge.
      int edge = hex_lattice_cell_edge(lattice, q, r, dir);
      mesh->cell_edges[mesh->cell_edge_offsets[cell]+dir] = edge;

      // Nodes for this edge.
      int node = hex_lattice_cell_node(lattice, q, r, dir);
      mesh->edge_nodes[2*edge] = node;
      mesh->edge_nodes[2*((edge+1)%6)+1] = node;

      // Cells to which this edge attaches.
      if (mesh->edge_cells[2*edge] == -1)
        mesh->edge_cells[2*edge] = cell;
      if (mesh->edge_cells[2*edge+1] == -1)
      {
        int q1, r1;
        hex_lattice_cell_get_neighbor(lattice, q, r, dir, &q1, &r1);
        int neighbor = hex_lattice_cell(lattice, q1, r1);
        mesh->edge_cells[2*edge] = neighbor;
      }
    }
  }

  // Clean up.
  hex_lattice_free(lattice);

  return mesh;
}

