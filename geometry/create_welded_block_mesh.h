// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_CREATE_WELDED_BLOCK_MESH_H
#define POLYMEC_CREATE_WELDED_BLOCK_MESH_H

#include "core/mesh.h"

// Given an array of logically-rectangular meshes (blocks) with boundary 
// (face) tags, this function welds together the meshes, identifying faces 
// with identical tags, to produce a single block-structured mesh containing 
// all of the cells, but only the nodes, edges, and faces that belong to the 
// intersection of the meshes. The nodal coordinates and mesh geometry of the 
// resulting mesh are the same as those of the input meshes. If two blocks 
// share a boundary tag but one of their faces or nodes is unmatched by another 
// within a distance given by the weld tolerance, this function fails fatally.
mesh_t* create_welded_block_mesh(mesh_t** blocks, int num_blocks, real_t weld_tolerance);

#endif

