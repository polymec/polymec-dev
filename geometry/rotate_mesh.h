// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_ROTATE_MESH_H
#define POLYMEC_ROTATE_MESH_H

#include "core/mesh.h"

// This function rotates the given mesh around the center point x0 by the 
// given angle theta about the axis defined by the pseudo-vector omega. 
// The node positions of the mesh are modified accordingly, but its geometry 
// is NOT recomputed.
void rotate_mesh(mesh_t* mesh, point_t* x0, real_t theta, vector_t* omega);

#endif

