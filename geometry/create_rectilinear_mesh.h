// Copyright 2012-2013 Jeffrey Johnson.
// 
// This file is part of Polymec, and is licensed under the Apache License, 
// Version 2.0 (the "License"); you may not use this file except in 
// compliance with the License. You may may find the text of the license in 
// the LICENSE file at the top-level source directory, or obtain a copy of 
// it at
// 
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef POLYMEC_CREATE_RECTILINEAR_MESH_H
#define POLYMEC_CREATE_RECTILINEAR_MESH_H

#include "core/mesh.h"

// This function creates and returns a rectilinear mesh whose nodes are 
// given by the xs, ys, and zs arrays.
mesh_t* create_rectilinear_mesh(MPI_Comm comm, 
                                double* xs, int nxs, 
                                double* ys, int nys, 
                                double* zs, int nzs);

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

