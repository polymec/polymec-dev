// Copyright (c) 2012-2017, Jeffrey N. Johnson
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
// ns - The number of cells on the surface-facing sides of each block.
// nr - The number of cells on the radial sides of each block.
// r - The radius of the (spherical) center block.
// R - The radius of the sphere being represented.
// R_tag - The name of the tag identifying faces along the surface of the sphere.
mesh_t* create_cubed_sphere_mesh(MPI_Comm comm,
                                 int ns, int nr, 
                                 real_t r, real_t R, 
                                 const char* R_tag);

// This function creates a mesh consisting of 6 blocks of hexahedral 
// cells in the shape of a spherical shell. The parameters of interest are:
// ns - The number of cells on the surface-facing sides of each block.
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

// This constructs a single panel of the cubed sphere with the given 
// properties. The panel is the top panel (centered at theta = 0) of the 
// sphere aligned with the z axis. The tags on this panel are prepended with 
// the given tag prefix, and are:
// - (tag_prefix)_0 - The faces that belong to the -x boundary.
// - (tag_prefix)_1 - The faces that belong to the +x boundary.
// - (tag_prefix)_2 - The faces that belong to the -y boundary.
// - (tag_prefix)_3 - The faces that belong to the +y boundary.
// - (tag_prefix)_4 - The faces that belong to the inner boundary.
// - (tag_prefix)_5 - The faces that belong to the outer boundary.
mesh_t* create_cubed_sphere_panel(MPI_Comm comm,
                                  int ns, int nr,
                                  real_t r, real_t R,
                                  const char* tag_prefix);

#endif

