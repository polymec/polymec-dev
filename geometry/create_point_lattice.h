// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_CREATE_POINT_LATTICE_H
#define POLYMEC_CREATE_POINT_LATTICE_H

#include "geometry/point_cloud.h"

/// \addtogroup geometry geometry
///@{

/// This function creates and returns a lattice of points whose x, y, and z
/// coordinates are given by arrays.
/// \relates point_cloud
point_cloud_t* create_point_lattice(MPI_Comm comm,
                                    real_t* xs, int nxs,
                                    real_t* ys, int nys,
                                    real_t* zs, int nzs);

/// This function creates and returns a uniform nx x ny x nz lattice of points.
/// The lattice spans the rectangular region of space defined by the bounding box.
/// \relates point_cloud
point_cloud_t* create_uniform_point_lattice(MPI_Comm comm, int nx, int ny, int nz, bbox_t* bbox);

///@}

#endif

