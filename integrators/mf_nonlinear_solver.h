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

#ifndef POLYMEC_MF_NONLINEAR_SOLVER_H
#define POLYMEC_MF_NONLINEAR_SOLVER_H

#include "kinsol/kinsol.h"
#include "core/adj_graph.h"
#include "integrators/nonlinear_solver.h"

// The different algorithms for matrix-free solution of nonlinear equations.
typedef enum
{
  GMRES,
  BICGSTAB,
  TFQMR
} mf_nonlinear_solver_type_t;

// Creates a solver that integrates the given set of nonlinear equations 
// using a preconditioned matrix-free Newton method.
nonlinear_solver_t* mf_nonlinear_solver_new(const char* name,
                                            void* context,
                                            KINSysFn F,
                                            nonlinear_solver_dtor dtor,
                                            adj_graph_t* graph,
                                            mf_nonlinear_solver_type_t type);

#endif

