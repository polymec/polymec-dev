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

#ifndef POLYMEC_POISSON_ELLIPTIC_SOLVER_H
#define POLYMEC_POISSON_ELLIPTIC_SOLVER_H

#include "core/st_func.h"
#include "core/mesh.h"
#include "integrators/elliptic_solver.h"

// Creates an elliptic solver for Poisson's equation with the given source 
// function, defined on the given mesh with the given set of boundary cells.
elliptic_solver_t* poisson_elliptic_solver_new(st_func_t* source,
                                               mesh_t* mesh,
                                               boundary_cell_map_t* boundary_cells,
                                               index_space_t* index_space);

#endif

