// Copyright (c) 2012-2014, Jeffrey N. Johnson
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this 
// list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice, 
// this list of conditions and the following disclaimer in the documentation 
// and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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
void tag_rectilinear_mesh_faces(mesh_t* mesh, int nx, int ny, int nz,
                                const char* x1_tag, 
                                const char* x2_tag, 
                                const char* y1_tag,
                                const char* y2_tag,
                                const char* z1_tag,
                                const char* z2_tag);

#endif

