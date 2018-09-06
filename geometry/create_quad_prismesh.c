// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "geometry/create_quad_planar_polymesh.h"
#include "geometry/create_quad_prismesh.h"

prismesh_t* create_quad_prismesh(MPI_Comm comm,
                                 size_t nx, size_t ny, size_t nz,
                                 bbox_t* bbox,
                                 bool periodic_in_x,
                                 bool periodic_in_y,
                                 bool periodic_in_z)
{
  planar_polymesh_t* columns = create_quad_planar_polymesh(nx, ny, bbox, 
                                                           periodic_in_x, 
                                                           periodic_in_y);
  return prismesh_new(comm, columns, bbox->z1, bbox->z2, nz, periodic_in_z);
}

