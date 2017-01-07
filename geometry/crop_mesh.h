// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_CROP_MESH_H
#define POLYMEC_CROP_MESH_H

#include "core/mesh.h"
#include "core/sp_func.h"

// Types of cropping algorithms.
typedef enum
{
  REMOVE_CELLS,  // cells are removed and nothing else is done
  PROJECT_NODES, // nodes of boundary cells are projected to the implicit function surface F(x) == 0.
  PROJECT_FACES  // boundary faces are placed in such a way to minimize the 2-norm of their distance
                 // from the boundary F(x) == 0.
} mesh_crop_t;

// This function marks the cells in a mesh whose centroids fall outside the 
// given implicit function, returning a copy of the mesh with these marked 
// cells removed. The resulting boundary faces of the cropped mesh are 
// tagged with the name of the boundary function.
mesh_t* crop_mesh(mesh_t* mesh, sp_func_t* boundary_func, mesh_crop_t crop_type);

#endif

