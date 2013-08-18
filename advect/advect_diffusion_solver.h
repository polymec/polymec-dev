// Copyright 2012-2013 Jeffrey Johnson.
// 
// This file is part of Polymec, and is licensed under the Apache License, 
// Version 2.0 (the "License"); you may not use this file except in 
// compliance with the License. You may may find the text of the license in 
// the LICENSE file at the top-level source directory, or obtain a copy of 
// it at
// 
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

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

