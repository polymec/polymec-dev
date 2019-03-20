// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_CREATE_RECTILINEAR_POLYMESH_H
#define POLYMEC_CREATE_RECTILINEAR_POLYMESH_H

#include "geometry/polymesh.h"

/// \addtogroup geometry geometry
///@{

/// This function creates and returns a rectilinear mesh whose nodes are
/// given by the xs, ys, and zs arrays. If comm == MPI_COMM_SELF, the
/// indices of the cells, faces, and nodes can all be navigated using a
/// an (nxs-1) x (nys-1) x (nzs-1) cubic lattice object.
/// \relates polymesh
polymesh_t* create_rectilinear_polymesh(MPI_Comm comm,
                                        real_t* xs, int nxs,
                                        real_t* ys, int nys,
                                        real_t* zs, int nzs);

/// This function creates and returns a rectilinear mesh whose nodes are
/// given by the xs, ys, and zs arrays on the given communicator
/// ONLY ON THE GIVEN RANK. The function returns the mesh on that rank and
/// a NULL pointer on all other ranks. MPI_COMM_SELF cannot be used as the
/// communicator.
/// \relates polymesh
polymesh_t* create_rectilinear_polymesh_on_rank(MPI_Comm comm,
                                                int rank,
                                                real_t* xs, int nxs,
                                                real_t* ys, int nys,
                                                real_t* zs, int nzs);

/// This function tags the faces of a rectilinear mesh for convenient boundary
/// condition assignments.
void tag_rectilinear_polymesh_faces(polymesh_t* mesh,
                                    const char* x1_tag,
                                    const char* x2_tag,
                                    const char* y1_tag,
                                    const char* y2_tag,
                                    const char* z1_tag,
                                    const char* z2_tag);

///@}

#endif

