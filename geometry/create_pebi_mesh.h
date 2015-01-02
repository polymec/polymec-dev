// Copyright (c) 2012-2015, Jeffrey N. Johnson
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

#ifndef POLYMEC_PEBI_H
#define POLYMEC_PEBI_H

#include "core/mesh.h"

// Perpendicular-bisector mesh features.
const char* PEBI;       // Is a PEBI mesh (has no edges or nodes)

// Creates a PEBI mesh given a set of cell centers and adjacency (face) 
// information. The faces array contains 2*num_faces entries, where 
// faces[2*i] and faces[2*i+1] contain the indices of the cells connected by 
// face i. If a face is only connected to one cell, faces[2*i+1] == -1.
// Meanwhile, face_areas contains num_faces entries, with face_areas[i] holding
// the area of face i. No edge or node information is stored.
// NOTE: If face_centers is NULL, each face center is assumed to lie at the 
// midpoint between its two cells. 
mesh_t* create_pebi_mesh(MPI_Comm comm, 
                         point_t* cell_centers, real_t* cell_volumes, int num_cells,
                         int* faces, real_t* face_areas, point_t* face_centers,
                         int num_faces);
 
// Creates a PEBI mesh from the given unstructured mesh. Features are 
// copied--tags are not.
mesh_t* create_pebi_mesh_from_unstructured_mesh(mesh_t* mesh);

#endif

