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
// given by the xs, ys, and zs arrays. If comm == MPI_COMM_SELF, the 
// indices of the cells, faces, and nodes can all be navigated using a
// an (nxs-1) x (nys-1) x (nzs-1) cubic lattice object.
mesh_t* create_rectilinear_mesh(MPI_Comm comm, 
                                real_t* xs, int nxs, 
                                real_t* ys, int nys, 
                                real_t* zs, int nzs);

// This function tags the faces of a rectilinear mesh for convenient boundary 
// condition assignments. This also creates a property named 
// "rectilinear_boundary_tags" containing an array of 6 strings:
// [x1_tag, x2_tag, y1_tag, y2_tag, z1_tag, z2_tag].
void tag_rectilinear_mesh_faces(mesh_t* mesh, 
                                const char* x1_tag, 
                                const char* x2_tag, 
                                const char* y1_tag,
                                const char* y2_tag,
                                const char* z1_tag,
                                const char* z2_tag);

#endif

