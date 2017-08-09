// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_MESH_STENCILS_H
#define POLYMEC_MESH_STENCILS_H

#include "core/mesh.h"
#include "model/stencil.h"

// Creates a star-shaped stencil for the cells in the given mesh. The stencil 
// has the given "radius," which is the maximum number of faces separating a 
// cell from one of its neighboring cells. This unweighted stencil is 
// constructed for every cell in the given mesh, and does not include the 
// "central" cell.
// NOTE: This stencil is not currently implemented for radius > 1!
stencil_t* cell_star_stencil_new(mesh_t* mesh, int radius);

// Creates a halo-shaped stencil for the cells in the given mesh. The halo 
// consists of a single layer of cells surrounding the nodes of each cell.
// This unweighted stencil is constructed for every cell in the given mesh, 
// and does not include the "central" cell.
// NOTE: This stencil is not currently implemented!
//stencil_t* cell_halo_stencil_new(mesh_t* mesh);

#endif
