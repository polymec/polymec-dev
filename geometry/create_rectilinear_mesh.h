// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_CREATE_RECTILINEAR_MESH_H
#define POLYMEC_CREATE_RECTILINEAR_MESH_H

#include "core/mesh.h"

// This function creates and returns a rectilinear mesh whose nodes are 
// given by the xs, ys, and zs arrays.
mesh_t* create_rectilinear_mesh(MPI_Comm comm, 
                                real_t* xs, int nxs, 
                                real_t* ys, int nys, 
                                real_t* zs, int nzs);

// This function tags the faces of a rectilinear mesh for convenient boundary 
// condition assignments.
void tag_rectilinear_mesh_faces(mesh_t* mesh, 
                                const char* x1_tag, 
                                const char* x2_tag, 
                                const char* y1_tag,
                                const char* y2_tag,
                                const char* z1_tag,
                                const char* z2_tag);

#endif

