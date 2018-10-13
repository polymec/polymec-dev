// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "geometry/create_hex_planar_polymesh.h"
#include "geometry/create_hex_prismesh.h"

prismesh_t* create_hex_prismesh(MPI_Comm comm,
                                size_t radius, real_t h, real_t rotation,
                                size_t nz, real_t z1, real_t z2,
                                bool periodic_in_z)
{
  ASSERT(h > 0.0);
  ASSERT(z1 < z2);
  ASSERT(nz > 0);
  planar_polymesh_t* columns = create_hex_planar_polymesh(radius, h, rotation);
  return prismesh_new(comm, columns, z1, z2, nz, periodic_in_z);
}

