// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_CREATE_VORONOI_MESH_H
#define POLYMEC_CREATE_VORONOI_MESH_H

#include "core/mesh.h"

// Creates a Voronoi mesh from a set of generator points within the given 
// bounding box using the algorithm described by Hugo Ledoux in his 2007 
// paper "Computing the 3D Voronoi Diagram Robustly: An Easy Explanation".
mesh_t* create_voronoi_mesh(MPI_Comm comm, point_t* generators, 
                            int num_generators, bbox_t* bounding_box);
 
#endif

