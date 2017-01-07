// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "geometry/translate_mesh.h"

void translate_mesh(mesh_t* mesh, vector_t* D)
{
  for (int n = 0; n < mesh->num_nodes; ++n)
  {
    point_t* xn = &mesh->nodes[n];
    xn->x += D->x;
    xn->y += D->y;
    xn->z += D->z;
  }
}

