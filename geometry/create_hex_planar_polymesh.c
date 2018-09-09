// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "geometry/create_hex_planar_polymesh.h"

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
  polymec_release(lattice);

  return mesh;
}

