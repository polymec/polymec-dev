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

#ifndef POLYMEC_ADVECT_DIFFUSION_SOLVER_H
#define POLYMEC_ADVECT_DIFFUSION_SOLVER_H

#include "core/st_func.h"
#include "core/mesh.h"
#include "integrators/diffusion_solver.h"

// Creates a diffusion solver model for the advection-diffusion equation 
// with the given diffusivity and source functions, defined on the given mesh
// with the given set of boundary cells.
diffusion_solver_t* advect_diffusion_solver_new(mesh_t* mesh,
                                                boundary_cell_map_t* boundary_cells,
                                                index_space_t* index_space);

// Resets the diffusivity function.
void advect_diffusion_solver_set_diffusivity(diffusion_solver_t* solver, st_func_t* diffusivity);

// Resets the source function.
void advect_diffusion_solver_set_source(diffusion_solver_t* solver, st_func_t* source);

// Sets the advective terms for use as source terms by the diffusion solver.
void advect_diffusion_solver_set_advective_source(diffusion_solver_t* solver,
                                                  double* advective_source);

#endif

