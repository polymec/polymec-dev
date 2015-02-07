// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_CREATE_CUBED_SPHERE_MESH_H
#define POLYMEC_CREATE_CUBED_SPHERE_MESH_H

#include "core/mesh.h"

// This function creates a mesh consisting of 7 blocks of hexahedral 
// cells in the shape of a cylinder. The parameters of interest are:
// n - The number of cells on the side of each of the 7 blocks.
// R - The radius of the sphere being represented.
// l - The length on a side of the center block of the sphere. This must 
//     be less than R.
// curved_center_block - If set to true, the center block is curved into a 
//                       sphere. If not, it is kept as a cube.
// R_tag - The name of the tag identifying faces along the surface of the sphere.
mesh_t* create_cubed_sphere_mesh(MPI_Comm comm,
                                 int n, real_t R, real_t l, 
                                 bool curved_center_block,
                                 const char* R_tag);

// This function creates a mesh consisting of 6 blocks of hexahedral 
// cells in the shape of a spherical shell. The parameters of interest are:
// na - The number of cells on the surface-facing sides of each block.
// nr - The number of cells on the radial sides of each block.
// r - The inner radius of the spherical shell.
// R - The outer radius of the spherical shell.
// The following tags are used to designate surfaces of interest:
// r_tag      - the set of faces along the inner radius of the shell.
// R_tag      - the set of faces along the outer radius of the shell.
mesh_t* create_cubed_spherical_shell_mesh(MPI_Comm comm,
                                          int ns, int nr,
                                          real_t r, real_t R,
                                          const char* r_tag,
                                          const char* R_tag);

// This type identifies a specific cubed sphere panel.
typedef enum
{
  X1_PANEL, // panel facing -x
  X2_PANEL, // panel facing +x
  Y1_PANEL, // panel facing -y
  Y2_PANEL, // panel facing +y
  Z1_PANEL, // panel facing -z
  Z2_PANEL  // panel facing +z
} cubed_sphere_panel_t;

// This constructs a single panel of the cubed sphere with the given properties,
// identified by which_panel.
mesh_t* create_cubed_sphere_panel(MPI_Comm comm,
                                  int ns, int nr,
                                  real_t R, real_t l,
                                  bool curved_bottom,
                                  cubed_sphere_panel_t which_panel);

#endif

