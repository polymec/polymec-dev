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

#ifndef POLYMEC_CREATE_DUAL_MESH_H
#define POLYMEC_CREATE_DUAL_MESH_H

#include "core/mesh.h"
#include "core/point.h"

// This function creates a mesh that is dual to the given original mesh, 
// respecting the features of the underlying geometric model, as identified 
// by the following tags: 
// - boundary_face_tags is an array containing names of face tags that denote 
//   interior or external boundaries within the model.
// - model_edge_tags is an array containing names of edge tags 
//   denoting edges conforming to the model. 
// - model_vertex_tags is an array containing names of node tags denoting 
//   vertices conforming to the model.
mesh_t* create_dual_mesh(MPI_Comm comm, 
                         mesh_t* original_mesh,
                         char** boundary_face_tags,
                         int num_boundary_face_tags,
                         char** model_edge_tags,
                         int num_model_edge_tags,
                         char** model_vertex_tags,
                         int num_model_vertex_tags);

#endif

