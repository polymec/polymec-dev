// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_CREATE_CUBED_CYLINDER_MESH_H
#define POLYMEC_CREATE_CUBED_CYLINDER_MESH_H

#include "core/mesh.h"

// This function creates a mesh consisting of 5 blocks of hexahedral 
// cells in the shape of a cylinder. The parameters of interest are:
// nx - The number of cells across a given block in the x-y plane.
// nz - The number of cells across a given block in the axial direction.
// R - The radius of the cylinder being represented.
// L - The length of the cylinder along the z axis.
// l - The length on a side of the middle block of the cylinder. This must 
//     be less than R.
// k - The curvature of the middle block, equal to 1/r, where r is the 
//     radius of the circle along which the middle block's boundary falls.
//     If k == 0, the middle block is an undeformed rectangular prism. If 
//     k == 1/(l*sqrt(2)), the middle block will itself be a cylinder with 
//     deformed cells.
mesh_t* create_cubed_cylinder_mesh(MPI_Comm comm,
                                   int nx, int nz,
                                   real_t R, real_t L,
                                   real_t l, real_t k);

// This function creates a mesh consisting of 4 blocks of hexahedral 
// cells in the shape of a cylindrical shell. The parameters of interest are:
// nx - The number of cells across a given block in the x-y plane.
// nz - The number of cells across a given block in the axial direction.
// r - The inner radius of the cylindrical shell.
// R - The outer radius of the cylindrical shell.
// L - The length of the cylinder along the z axis.
mesh_t* create_cubed_cylindrical_shell_mesh(MPI_Comm comm,
                                            int nx, int nz,
                                            real_t r, real_t R, real_t L);

#endif

