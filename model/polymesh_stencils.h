// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_POLYMESH_STENCILS_H
#define POLYMEC_POLYMESH_STENCILS_H

#include "geometry/polymesh.h"
#include "model/stencil.h"

/// \addtogroup model model
///@{

/// Creates a star-shaped stencil for the cells in the given mesh. The stencil
/// has the given "radius," which is the maximum number of faces separating a
/// cell from one of its neighboring cells. This unweighted stencil is
/// constructed for every cell in the given mesh, and does not include the
/// "central" cell.
/// \note This stencil is not currently implemented for radius > 1!
stencil_t* cell_star_stencil_new(polymesh_t* mesh, int radius);

///@}

#endif
