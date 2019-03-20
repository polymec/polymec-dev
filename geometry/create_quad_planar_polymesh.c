// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "geometry/create_quad_planar_polymesh.h"

typedef struct
{
  // Number of cells in x, y, and z.
  size_t nx, ny;
} quad_lattice_t;

static quad_lattice_t* quad_lattice_new(size_t nx, size_t ny)
{
  ASSERT(nx > 0);
  ASSERT(ny > 0);
  quad_lattice_t* l = polymec_malloc(sizeof(quad_lattice_t));
  l->nx = nx;
  l->ny = ny;
  return l;
}

static void quad_lattice_free(quad_lattice_t* l)
{
  polymec_free(l);
}

static inline size_t quad_lattice_num_cells(quad_lattice_t* l)
{
  return l->nx * l->ny;
}

static inline size_t quad_lattice_num_x_edges(quad_lattice_t* l)
{
  return l->nx * (l->ny+1);
}

static inline size_t quad_lattice_num_y_edges(quad_lattice_t* l)
{
  return (l->nx+1) * l->ny;
}

static inline size_t quad_lattice_num_edges(quad_lattice_t* l)
{
  return quad_lattice_num_x_edges(l) + quad_lattice_num_y_edges(l);
}

static inline size_t quad_lattice_num_nodes(quad_lattice_t* l)
{
  return (l->nx+1) * (l->ny+1);
}

static inline int quad_lattice_cell(quad_lattice_t* l, int i, int j)
{
  return (int)(l->nx*j + i);
}

static inline void quad_lattice_get_cell_pair(quad_lattice_t* l, int index, int* i, int* j)
{
  *j = (int)(index/l->nx);
  *i = (int)(index - l->nx*(*j));
}

static inline int quad_lattice_x_edge(quad_lattice_t* l, int i, int j)
{
  return (int)l->nx*j + i;
}

static inline int quad_lattice_y_edge(quad_lattice_t* l, int i, int j)
{
  return (int)(quad_lattice_num_x_edges(l) + (l->nx+1)*j + i);
}

static inline int quad_lattice_node(quad_lattice_t* l, int i, int j)
{
  return (int)((l->nx+1)*j + i);
}

planar_polymesh_t* create_quad_planar_polymesh(size_t nx, size_t ny,
                                               bbox_t* bbox,
                                               bool periodic_in_x,
                                               bool periodic_in_y)
{
  ASSERT(bbox->x2 > bbox->x1);
  ASSERT(bbox->y2 > bbox->y1);

  // Create a cubic lattice object for indexing.
  quad_lattice_t* lattice = quad_lattice_new(nx, ny);

  int num_cells = (int)quad_lattice_num_cells(lattice);
  int num_edges = (int)quad_lattice_num_edges(lattice);
  int num_nodes = (int)quad_lattice_num_nodes(lattice);
  int num_edges_per_cell = 4;
  planar_polymesh_t* mesh = planar_polymesh_new_with_cell_type(num_cells,
                                                               num_edges,
                                                               num_nodes,
                                                               num_edges_per_cell);

  // Now traverse the cells and use the lattice to assign edges and nodes.
  for (int c = 0; c < num_cells; ++c)
  {
    // Figure out (i, j) indices for this cell.
    int cell_index = c, i, j;
    quad_lattice_get_cell_pair(lattice, cell_index, &i, &j);

    // Hook up the cell's edges.
    int edges[4] = {quad_lattice_x_edge(lattice, i,     j),
                    quad_lattice_y_edge(lattice, i+1,   j),
                    quad_lattice_x_edge(lattice, i,   j+1),
                    quad_lattice_y_edge(lattice, i,     j)};
    for (int e = 0; e < 4; ++e)
    {
      if (mesh->edge_cells[2*edges[e]] == -1)
      {
        mesh->cell_edges[mesh->cell_edge_offsets[c] + e] = edges[e];
        mesh->edge_cells[2*edges[e]] = c;
      }
      else
      {
        // An edge is attached to at most 2 cells.
        ASSERT(mesh->edge_cells[2*edges[e]+1] == -1);
        mesh->cell_edges[mesh->cell_edge_offsets[c] + e] = ~(edges[e]);
        mesh->edge_cells[2*edges[e]+1] = c;
      }
    }

    // Hook up each edge's nodes.
    int nodes[4] = {quad_lattice_node(lattice, i,   j),
                    quad_lattice_node(lattice, i+1, j),
                    quad_lattice_node(lattice, i+1, j+1),
                    quad_lattice_node(lattice, i,   j+1)};
    for (int n = 0; n < 4; ++n)
    {
      mesh->edge_nodes[2*edges[n]] = nodes[n];
      mesh->edge_nodes[2*edges[n]+1] = nodes[(n+1)%4];
    }
  }

  // Set node positions.
  real_t dx = (bbox->x2 - bbox->x1) / nx;
  real_t dy = (bbox->y2 - bbox->y1) / ny;
  for (int i = 0; i <= (int)nx; ++i)
  {
    for (int j = 0; j <= (int)ny; ++j)
    {
      int n = quad_lattice_node(lattice, i, j);
      mesh->nodes[n].x = i*dx;
      mesh->nodes[n].y = j*dy;
    }
  }

  // Clean up.
  quad_lattice_free(lattice);

  return mesh;
}

