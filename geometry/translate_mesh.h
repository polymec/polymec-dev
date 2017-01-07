// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_TRANSLATE_MESH_H
#define POLYMEC_TRANSLATE_MESH_H

#include "core/mesh.h"

// This function displaces all of the nodes of the given mesh by the given 
// displacement vector D. The geometry is not recomputed, so face and cell 
// centers will be incorrect until it is.
void translate_mesh(mesh_t* mesh, vector_t* D);

#endif

