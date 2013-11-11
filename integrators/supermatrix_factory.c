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

#include "integrators/supermatrix_factory.h"

struct supermatrix_factory_t 
{
  adj_graph_t* graph;
  adj_graph_coloring_t* coloring;
  KINSysFn F;
  CVRhsFn rhs;
  void* context;
};

supermatrix_factory_t* supermatrix_factory_from_sys_func(adj_graph_t* graph,
                                                         KINSysFn F,
                                                         void* context)
{
  supermatrix_factory_t* factory = malloc(sizeof(supermatrix_factory_t));
  factory->graph = graph;
  factory->coloring = adj_graph_coloring_new(graph, SMALLEST_LAST);
  factory->F = F;
  factory->rhs = NULL;
  factory->context = context;
  return factory;
}

supermatrix_factory_t* supermatrix_factory_from_rhs(adj_graph_t* graph,
                                                    CVRhsFn rhs,
                                                    void* context)
{
  supermatrix_factory_t* factory = malloc(sizeof(supermatrix_factory_t));
  factory->graph = graph;
  factory->coloring = adj_graph_coloring_new(graph, SMALLEST_LAST);
  factory->F = NULL;
  factory->rhs = rhs;
  factory->context = context;
  return factory;
}

void supermatrix_factory_free(supermatrix_factory_t* factory)
{
  adj_graph_coloring_free(factory->coloring);
  free(factory);
}

SuperMatrix* supermatrix_factory_jacobian(supermatrix_factory_t* factory, N_Vector u, double t)
{
  return NULL;
}

void supermatrix_factory_update_jacobian(supermatrix_factory_t* factory, N_Vector u, double t, SuperMatrix* J)
{
}

