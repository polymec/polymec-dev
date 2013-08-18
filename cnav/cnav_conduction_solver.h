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

#ifndef POLYMEC_CNAV_CONDUCTION_SOLVER_H
#define POLYMEC_CNAV_CONDUCTION_SOLVER_H

#include "core/diffusion_solver.h"
#include "core/st_func.h"
#include "core/mesh.h"

// Creates a diffusion solver model for heat conduction in the 
// compressible Navier-Stokes model.
diffusion_solver_t* cnav_conduction_solver_new(mesh_t* mesh,
                                               boundary_cell_map_t* boundary_cells);

// Sets the cnavive terms for use as source terms by the diffusion solver.
void cnav_conduction_solver_set_advective_source(diffusion_solver_t* solver,
                                                 double* advective_source);

#endif

