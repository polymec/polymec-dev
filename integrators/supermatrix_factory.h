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

#ifndef POLYMEC_SUPERMATRIX_FACTORY_H
#define POLYMEC_SUPERMATRIX_FACTORY_H

#include "core/polymec.h"
#include "core/adj_graph.h"
#include "cvode/cvode.h"
#include "kinsol/kinsol.h"
#include "slu_ddefs.h"
#include "supermatrix.h"

// This class represents a way of integrating a system of 
// nonlinear equations.
typedef struct supermatrix_factory_t supermatrix_factory_t;

// Creates a factory that produces SuperLU SuperMatrix objects using 
// the given adjacency graph and system function F (and context). The 
// factory does NOT assume ownership of the factory or context.
supermatrix_factory_t* supermatrix_factory_from_sys_func(adj_graph_t* graph,
                                                         KINSysFn F,
                                                         void* context);

// Creates a factory that produces SuperLU SuperMatrix objects using 
// the given adjacency graph and RHS function (and context). The
// factory does NOT assume ownership of the factory or context.
supermatrix_factory_t* supermatrix_factory_from_rhs(adj_graph_t* graph,
                                                    CVRhsFn rhs,
                                                    void* context);

// Frees the factory.
void supermatrix_factory_free(supermatrix_factory_t* factory);

// Produces a supermatrix with a sparsity that matches the given graph 
// from which the factory was constructed.
SuperMatrix* supermatrix_factory_matrix(supermatrix_factory_t* factory);

// Produces a Jacobian supermatrix from finite difference quotients applied 
// to the given function F at the solution u and t.
SuperMatrix* supermatrix_factory_jacobian(supermatrix_factory_t* factory, N_Vector u, double t);

// Updates a Jacobian supermatrix from finite difference quotients applied 
// to the given function F at the solution u and t.
void supermatrix_factory_update_jacobian(supermatrix_factory_t* factory, N_Vector u, double t, SuperMatrix* J);

// Call this function to destroy supermatrices that have been created by this 
// factory.
void supermatrix_free(SuperMatrix* matrix);

#endif

