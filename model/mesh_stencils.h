// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_MESH_STENCILS_H
#define POLYMEC_MESH_STENCILS_H

#include "core/mesh.h"
#include "model/stencil.h"

// Creates a stencil for the cells that share at least one face with a given 
// cell. The stencil is constructed for every cell in the given mesh.
stencil_t* cell_face_stencil_new(mesh_t* mesh);

// Creates a stencil for the cells that share at least one edge with a given 
// cell. The stencil is constructed for every cell in the given mesh.
stencil_t* cell_edge_stencil_new(mesh_t* mesh);

// Creates a stencil for the cells that share at least one node with a given 
// cell. The stencil is constructed for every cell in the given mesh.
stencil_t* cell_node_stencil_new(mesh_t* mesh);

#endif
