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
// l - The length on a side of the center block of the cylinder. This must 
//     be less than R.
// curved_center_block - If set to true, the center block is curved into a 
//                       circle. If not, it is kept as a square.
// The following tags are used to designate surfaces of interest:
// R_tag      - the set of faces along the radius of the cylinder.
// bottom_tag - the set of faces along the bottom of the cylinder.
// top_tag    - the set of faces along the top of the cylinder.
mesh_t* create_cubed_cylinder_mesh(MPI_Comm comm,
                                   int nx, int nz,
                                   real_t R, real_t L, real_t l, 
                                   bool curved_center_block,
                                   const char* R_tag,
                                   const char* bottom_tag,
                                   const char* top_tag);

// This function creates a mesh consisting of 4 blocks of hexahedral 
// cells in the shape of a cylindrical shell. The parameters of interest are:
// nx - The number of cells across a given block in the x-y plane.
// nz - The number of cells across a given block in the axial direction.
// r - The inner radius of the cylindrical shell.
// R - The outer radius of the cylindrical shell.
// L - The length of the cylinder along the z axis.
// The following tags are used to designate surfaces of interest:
// r_tag      - the set of faces along the inner radius of the shell.
// R_tag      - the set of faces along the outer radius of the shell.
// bottom_tag - the set of faces along the bottom of the shell.
// top_tag    - the set of faces along the top of the shell.
mesh_t* create_cubed_cylindrical_shell_mesh(MPI_Comm comm,
                                            int nx, int nz,
                                            real_t r, real_t R, real_t L,
                                            const char* r_tag,
                                            const char* R_tag,
                                            const char* bottom_tag,
                                            const char* top_tag);

#endif

