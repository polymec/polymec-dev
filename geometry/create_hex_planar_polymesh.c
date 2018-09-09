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

  // Start at (0, 0) and spiral outward.
  // FIXME

  // Clean up.
  polymec_release(lattice);

  return mesh;
}

