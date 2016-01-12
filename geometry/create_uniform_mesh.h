// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_CREATE_UNIFORM_MESH_H
#define POLYMEC_CREATE_UNIFORM_MESH_H

#include "core/point.h"
#include "core/mesh.h"
#include "geometry/create_rectilinear_mesh.h"

// This function creates and returns a uniform mesh of nx x ny x nz cells. 
// The mesh spans the rectangular region of space defined by the bounding box.
mesh_t* create_uniform_mesh(MPI_Comm comm, int nx, int ny, int nz, bbox_t* bbox);

// This function creates and returns a uniform mesh of nx x ny x nz cells
// on the given communicator ONLY ON THE GIVEN RANK. The function returns the 
// mesh on that rank and a NULL pointer on all other ranks. MPI_COMM_SELF 
// cannot be used as the communicator.
mesh_t* create_uniform_mesh_on_rank(MPI_Comm comm, int rank,
                                    int nx, int ny, int nz, bbox_t* bbox);
#endif

