// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_CREATE_UNIFORM_POLYMESH_H
#define POLYMEC_CREATE_UNIFORM_POLYMESH_H

#include "geometry/create_rectilinear_polymesh.h"

/// \addtogroup geometry geometry
///@{

/// This function creates and returns a uniform mesh of nx x ny x nz cells.
/// The mesh spans the rectangular region of space defined by the bounding box.
/// \relates polymesh
polymesh_t* create_uniform_polymesh(MPI_Comm comm, int nx, int ny, int nz, bbox_t* bbox);

/// This function creates and returns a uniform mesh of nx x ny x nz cells
/// on the given communicator ONLY ON THE GIVEN RANK. The function returns the
/// mesh on that rank and a NULL pointer on all other ranks. MPI_COMM_SELF
/// cannot be used as the communicator, but it will ultimately be the
/// communicator associated with the mesh returned.
/// \relates polymesh
polymesh_t* create_uniform_polymesh_on_rank(MPI_Comm comm, int rank,
                                            int nx, int ny, int nz, bbox_t* bbox);

///@}

#endif

