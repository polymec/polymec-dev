// Copyright (c) 2012-2017, Jeffrey N. Johnson
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
    vector_t Y1;
    Y1.x = Y.x * cos_theta - Y.y * sin_theta;
    Y1.y = Y.x * sin_theta + Y.y * cos_theta;
    Y1.z = Y.z;

    // Project back.
    vector_t X1;
    X1.x = Y1.x * e1.x + Y1.y * e2.x + Y1.z * e3.x;
    X1.y = Y1.x * e1.y + Y1.y * e2.y + Y1.z * e3.y;
    X1.z = Y1.x * e1.z + Y1.y * e2.z + Y1.z * e3.z;

    mesh->nodes[n].x = x0->x + X1.x;
    mesh->nodes[n].y = x0->y + X1.y;
    mesh->nodes[n].z = x0->z + X1.z;
  }
}

