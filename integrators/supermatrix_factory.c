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
#include "slu_util.h"

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

SuperMatrix* supermatrix_factory_matrix(supermatrix_factory_t* factory)
{
  SuperMatrix* A = malloc(sizeof(SuperMatrix));
  
  // Fetch sparsity information from the graph.
  adj_graph_t* graph = factory->graph;
  int* edge_offsets = adj_graph_edge_offsets(graph);
  int num_rows = adj_graph_num_vertices(graph);
  int num_edges = edge_offsets[num_rows+1];
  int num_nz = num_rows + num_edges;

  // Translate the sparsity information in the graph to column indices and 
  // row pointers. Since the graph is effectively stored in compressed 
  // row format, we will use SuperLU's compressed row matrix format.
  int* row_ptrs = intMalloc(num_rows + 1);
  memcpy(row_ptrs, edge_offsets, sizeof(int) * (num_rows+1));
  int* col_indices = intMalloc(num_nz);
  int* edges = adj_graph_adjacency(graph);
  int offset = 0;
  for (int i = 0; i < num_rows; ++i)
  {
    col_indices[offset++] = i; // diagonal column
    for (int j = edge_offsets[i]; j < edge_offsets[i+1]; ++j)
      col_indices[offset++] = edges[j];
  }

  // Create zeros for the matrix.
  double* mat_zeros = doubleMalloc(num_nz);
  memset(mat_zeros, 0, sizeof(double) * num_nz);

  // Hand over these resources to create the Supermatrix.
  dCreate_CompRow_Matrix(A, num_rows, num_rows, num_nz, 
                         mat_zeros, col_indices, row_ptrs, 
                         SLU_NR, SLU_D, SLU_GE);

  return A;
}

SuperMatrix* supermatrix_factory_jacobian(supermatrix_factory_t* factory, N_Vector u, double t)
{
  SuperMatrix* J = supermatrix_factory_matrix(factory);
  supermatrix_factory_update_jacobian(factory, u, t, J);
  return J;
}

void supermatrix_factory_update_jacobian(supermatrix_factory_t* factory, N_Vector u, double t, SuperMatrix* J)
{
  // FIXME
}

void supermatrix_free(SuperMatrix* matrix)
{
  Destroy_CompCol_Matrix(matrix);
  free(matrix);
}

