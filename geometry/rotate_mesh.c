// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "geometry/rotate_mesh.h"

void rotate_mesh(mesh_t* mesh, point_t* x0, real_t theta, vector_t* omega)
{
  // Form a 3D coordinate basis given the axis.
  vector_t e1, e2, e3 = *omega;
  compute_orthonormal_basis(&e3, &e1, &e2); // e1 x e2 == e3.

  // Now rotate the nodes about x0.
  real_t cos_theta = cos(theta);
  real_t sin_theta = sin(theta);
  for (int n = 0; n < mesh->num_nodes; ++n)
  {
    // Express the node position in our new basis, centered about x0.
    vector_t X, Y;
    point_displacement(x0, &mesh->nodes[n], &X);
    Y.x = vector_dot(&X, &e1);
    Y.y = vector_dot(&X, &e2);
    Y.z = vector_dot(&X, &e3);

    // Rotate.
    Y.x =  Y.x * cos_theta + Y.y * sin_theta;
    Y.y = -Y.x * sin_theta + Y.y * cos_theta;

    // Project back.
    X.x = Y.x * e1.x + Y.y * e2.x + Y.z * e3.x;
    X.y = Y.x * e1.y + Y.y * e2.y + Y.z * e3.y;
    X.z = Y.x * e1.z + Y.y * e2.z + Y.z * e3.z;

    mesh->nodes[n].x = x0->x + X.x;
    mesh->nodes[n].y = x0->y + X.y;
    mesh->nodes[n].z = x0->z + X.z;
  }
}

