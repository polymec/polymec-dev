// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_CREATE_DUAL_MESH_H
#define POLYMEC_CREATE_DUAL_MESH_H

#include "core/mesh.h"
#include "core/point.h"

// This function creates a mesh that is dual to the given original mesh, 
// respecting the features of the underlying geometric model, as identified 
// by the following tags: 
// - external_model_face_tags is an array containing names of face tags that 
//   denote external boundaries within the model.
// - internal_model_face_tags is an array containing names of face tags that 
//   denote internal boundaries (separating, say, different materials) 
//   within the model.
// - model_edge_tags is an array containing names of edge tags 
//   denoting edges conforming to the model. 
// - model_vertex_tags is an array containing names of node tags denoting 
//   vertices conforming to the model.
mesh_t* create_dual_mesh(MPI_Comm comm, 
                         mesh_t* original_mesh,
                         char** external_model_face_tags,
                         int num_external_model_face_tags,
                         char** internal_model_face_tags,
                         int num_internal_model_face_tags,
                         char** model_edge_tags,
                         int num_model_edge_tags,
                         char** model_vertex_tags,
                         int num_model_vertex_tags);

#endif

