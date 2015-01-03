// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_CREATE_TETGEN_MESH_H
#define POLYMEC_CREATE_TETGEN_MESH_H

#include "core/mesh.h"
#include "core/point.h"

// This function creates a mesh using .node, .ele, .face, and .neigh files 
// created by TetGen. A global mesh is constructed on process 0 and is then 
// partitioned onto the given communicator using an unweighted partitioning 
// algorithm with a maximum load imbalance ratio of 0.05.
mesh_t* create_tetgen_mesh(MPI_Comm comm, 
                           const char* node_file,
                           const char* ele_file,
                           const char* face_file,
                           const char* neigh_file);

#endif

