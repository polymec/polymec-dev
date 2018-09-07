// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "geometry/create_quad_planar_polymesh.h"
#include "geometry/quad_lattice.h"

planar_polymesh_t* create_quad_planar_polymesh(size_t nx, size_t ny, 
                                               bbox_t* bbox,
                                               bool periodic_in_x,
                                               bool periodic_in_y)
{
  ASSERT(bbox->x2 > bbox->x1);
  ASSERT(bbox->y2 > bbox->y1);
  
  int num_cells = (int)(nx * ny);
  int num_edges = (int)((nx+1)*ny + nx*(ny+1));
  int num_nodes = (int)((nx+1) * (ny+1));
  int num_edges_per_cell = 4;
  planar_polymesh_t* mesh = planar_polymesh_new_with_cell_type(num_cells,
                                                               num_edges,
                                                               num_nodes,
                                                               num_edges_per_cell);

  // Create a cubic lattice object for indexing.
  quad_lattice_t* lattice = quad_lattice_new(nx, ny);

  // Now traverse the cells and use the lattice to assign edges and nodes.
  for (int c = 0; c < num_cells; ++c)
  {
    // Figure out (i, j) indices for this cell.
    index_t cell_index = (index_t)c, i, j;
    quad_lattice_get_cell_pair(lattice, cell_index, &i, &j);

    // Edges are just xy faces, really. Count em up.
    index_t edges[4] = {quad_lattice_x_face(lattice, i, j),
                        quad_lattice_x_face(lattice, i+1, j),
                        quad_lattice_y_face(lattice, i, j),
                        quad_lattice_y_face(lattice, i, j+1)};
    for (int e = 0; e < 4; ++e)
    {
      mesh->cell_edges[mesh->cell_edge_offsets[c] + e] = (int)edges[e];
      if (mesh->edge_cells[2*edges[e]] == -1)
        mesh->edge_cells[2*edges[e]] = c;
      else
        mesh->edge_cells[2*edges[e]+1] = c;
    }

    // We use the "bottom" nodes for the cell.
    index_t nodes[4] = {quad_lattice_node(lattice, i,   j),
                        quad_lattice_node(lattice, i+1, j),
                        quad_lattice_node(lattice, i+1, j+1),
                        quad_lattice_node(lattice, i,   j+1)};
    for (int n = 0; n < 4; ++n)
    {
      mesh->edge_nodes[2*edges[n]] = (int)nodes[n];
      mesh->edge_nodes[2*edges[(n+1)%4]] = (int)nodes[(n+1)%4];
    }
  }

  // Set node positions.
  real_t dx = (bbox->x2 - bbox->x1) / nx;
  real_t dy = (bbox->y2 - bbox->y1) / ny;
  for (index_t i = 0; i <= (index_t)nx; ++i)
  {
    for (index_t j = 0; j <= (index_t)ny; ++j)
    {
      index_t n = quad_lattice_node(lattice, i, j);
      mesh->nodes[n].x = i*dx;
      mesh->nodes[n].y = j*dy;
    }
  }

  return mesh;
}

