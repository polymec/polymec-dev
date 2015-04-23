// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_SYMMETRY_H
#define POLYMEC_SYMMETRY_H

#include "core/mesh.h"

// Symmetry-related mesh features.
extern const char* SYMMETRIC;       // Has 1D or 2D symmetry.
extern const char* ONE_DIMENSIONAL; // Has 1D symmetry.
extern const char* TWO_DIMENSIONAL; // Has 2D symmetry.
extern const char* CARTESIAN_1D;    // Has 1D linear Cartesian geometry.
extern const char* CARTESIAN_2D;    // Has 2D planar Cartesian geometry.
extern const char* CYLINDRICAL_1D;  // Has 1D radial cylindrical symmetry.
extern const char* CYLINDRICAL_RZ;  // Has 2D r-z cylindrical symmetry.
extern const char* SPHERICAL_1D;    // Has 1D radial spherical symmetry.

// These functions create meshes that can be used to run 1D problems.

// Creates a uniform line of N Cartesian cells spanning a domain [x1, x2].
mesh_t* create_uniform_cartesian_1d_mesh(MPI_Comm comm, real_t x1, real_t x2, int N);

// Creates a line of N Cartesian cells spanning [x1, x2], with cell lengths that 
// grow according to the given logarithmic factor. The length of the (i+1)th 
// cell will be log_factor times that of the ith cell.
mesh_t* create_logarithmic_cartesian_1d_mesh(MPI_Comm comm, real_t x1, real_t x2, real_t log_factor, int N);

// Creates a nonuniform line of N Cartesian cells with the given nodal 
// spacings xs. The array xs should be N+1 in length.
mesh_t* create_nonuniform_cartesian_1d_mesh(MPI_Comm comm, real_t* xs, int N);

// Creates a uniform line of N cylindrical cells spanning a domain [r1, r2].
mesh_t* create_uniform_cylindrical_1d_mesh(MPI_Comm comm, real_t r1, real_t r2, int N);

// Creates a line of N cylindrical cells spanning [r1, r2], with cell lengths 
// that grow according to the given logarithmic factor. The length of the 
// (i+1)th cell will be log_factor times that of the ith cell.
mesh_t* create_logarithmic_cylindrical_1d_mesh(MPI_Comm comm, real_t r1, real_t r2, real_t log_factor, int N);

// Creates a nonuniform line of N cylindrical cells with the given nodal 
// spacings rs. The array rs should be N+1 in length.
mesh_t* create_nonuniform_cylindrical_1d_mesh(MPI_Comm comm, real_t* rs, int N);

// Creates a uniform line of N spherical cells spanning a domain [r1, r2].
mesh_t* create_uniform_spherical_1d_mesh(MPI_Comm comm, real_t r1, real_t r2, int N);

// Creates a line of N spherical cells spanning [r1, r2] with cell lengths 
// that grow according to the given logarithmic factor. The length of the 
// (i+1)th cell will be log_factor times that of the ith cell.
mesh_t* create_logarithmic_spherical_1d_mesh(MPI_Comm comm, real_t r1, real_t r2, real_t log_factor, int N);

// Creates a nonuniform line of N spherical cells with the given nodal 
// spacings rs. The array rs should be N+1 in length.
mesh_t* create_nonuniform_spherical_1d_mesh(MPI_Comm comm, real_t* rs, int N);

#endif

