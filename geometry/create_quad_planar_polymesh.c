// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "geometry/create_quad_planar_polymesh.h"

planar_polymesh_t* create_quad_planar_polymesh(int nx, int ny, 
                                               bbox_t* bbox,
                                               bool periodic_in_x,
                                               bool periodic_in_y)
{
  int num_cells = nx * ny;
  int num_edges = (nx+1)*ny + nx*(ny+1);
  int num_nodes = (nx+1) * (ny+1);
  int num_edges_per_cell = 4;
  planar_polymesh_t* mesh = planar_polymesh_new_with_cell_type(num_cells,
                                                               num_edges,
                                                               num_nodes,
                                                               num_edges_per_cell);
  return mesh;
}

