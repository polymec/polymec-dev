// Copyright (c) 2012-2013, Jeffrey N. Johnson
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this 
// list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice, 
// this list of conditions and the following disclaimer in the documentation 
// and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef POLYMEC_SYMMETRY_1D_H
#define POLYMEC_SYMMETRY_1D_H

#include "core/mesh.h"

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

