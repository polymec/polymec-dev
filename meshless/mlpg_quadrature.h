// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef MLPG_QUADRATURE_H
#define MLPG_QUADRATURE_H

#include "core/point_cloud.h"
#include "integrators/volume_integral.h"
#include "integrators/surface_integral.h"

// This file contains several quadrature rules for evaluating volume and 
// surface integrals using the Meshless Local Petrov-Galerkin method, which 
// does not rely on a background mesh.

// This quadrature rule computes surface integrals over cubes whose extents 
// are defined by the ratio of its side to the extent associated with a 
// subdomain. The rule is associated with the given point cloud and (scalar) 
// subdomain extents field, and has the given number of points on a side.
surface_integral_t* mlpg_cube_surface_integral_new(point_cloud_t* cloud,
                                                   real_t* extents,
                                                   int num_points,
                                                   real_t side_to_extent_ratio);

// This quadrature rule computes volume integrals over cubes whose extents 
// are defined by the ratio of its side to the extent associated with a 
// subdomain. The rule is associated with the given point cloud and (scalar) 
// subdomain extents field, and has the given number of points on a side.
volume_integral_t* mlpg_cube_volume_integral_new(point_cloud_t* cloud,
                                                 real_t* extents,
                                                 int num_points,
                                                 real_t side_to_extent_ratio);

// This quadrature rule computes surface integrals over spheres whose extents 
// are defined by the ratio of its radius to the extent associated with a 
// subdomain. The rule is associated with the given point cloud and (scalar) 
// subdomain extents field, and has the given number of azimuthal points.
surface_integral_t* mlpg_sphere_surface_integral_new(point_cloud_t* cloud,
                                                     real_t* extents,
                                                     int num_points,
                                                     real_t radius_to_extent_ratio);

// This quadrature rule computes volume integrals over spheres whose extents 
// are defined by the ratio of its radius to the extent associated with a 
// subdomain. The rule is associated with the given point cloud and (scalar) 
// subdomain extents field, and has the given number of radial points.
volume_integral_t* mlpg_sphere_volume_integral_new(point_cloud_t* cloud,
                                                   real_t* extents,
                                                   int num_points,
                                                   real_t radius_to_extent_ratio);

#endif
