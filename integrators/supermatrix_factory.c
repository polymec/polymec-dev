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

#include "core/array_utils.h"
#include "core/sundials_helpers.h"
#include "integrators/supermatrix_factory.h"
#include "slu_util.h"

struct supermatrix_factory_t 
{
  adj_graph_t* graph;
  adj_graph_coloring_t* coloring;
  KINSysFn F;
  void (*set_F_time)(double, void*);
  CVRhsFn rhs;
  void* context;

  // Work vectors.
  N_Vector* work;
};

supermatrix_factory_t* supermatrix_factory_from_sys_func(adj_graph_t* graph,
                                                         KINSysFn F,
                                                         void (*set_F_time)(double, void*),
                                                         void* context)
{
  supermatrix_factory_t* factory = malloc(sizeof(supermatrix_factory_t));
  factory->graph = graph;
  factory->coloring = adj_graph_coloring_new(graph, SMALLEST_LAST);
  factory->F = F;
  factory->set_F_time = set_F_time;
  factory->rhs = NULL;
  factory->context = context;

  // Make work vectors.
  MPI_Comm comm = adj_graph_comm(factory->graph);
  int num_vertices = adj_graph_num_vertices(factory->graph);
  N_Vector prototype = N_VNew(comm, num_vertices);
//  int num_colors = adj_graph_coloring_num_colors(factory->coloring);
//  factory->work = N_VCloneVectorArray(num_colors, prototype);
  factory->work = N_VCloneVectorArray(4, prototype);
  N_VDestroy(prototype);
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
  factory->set_F_time = NULL;
  factory->rhs = rhs;
  factory->context = context;

  // Make work vectors.
  MPI_Comm comm = adj_graph_comm(factory->graph);
  int num_vertices = adj_graph_num_vertices(factory->graph);
  N_Vector prototype = N_VNew(comm, num_vertices);
//  int num_colors = adj_graph_coloring_num_colors(factory->coloring);
//  factory->work = N_VCloneVectorArray(num_colors, prototype);
  factory->work = N_VCloneVectorArray(4, prototype);
  N_VDestroy(prototype);
  return factory;
}

void supermatrix_factory_free(supermatrix_factory_t* factory)
{
//  int num_colors = adj_graph_coloring_num_colors(factory->coloring);
//  N_VDestroyVectorArray(factory->work, num_colors);
  N_VDestroyVectorArray(factory->work, 4);
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
  int num_edges = edge_offsets[num_rows];
  int num_nz = num_rows + num_edges;

  // Translate the sparsity information in the graph to column indices and 
  // row pointers. Since the graph is effectively stored in compressed 
  // row format, we will use SuperLU's compressed row matrix format.
  int* row_ptrs = intMalloc(num_rows + 1);
  int* col_indices = intMalloc(num_nz);
  int* edges = adj_graph_adjacency(graph);
  int offset = 0;
  for (int i = 0; i < num_rows; ++i)
  {
    row_ptrs[i] = offset;
    col_indices[offset++] = i; // diagonal column

    // Off-diagonal columns.
    for (int j = edge_offsets[i]; j < edge_offsets[i+1]; ++j)
      col_indices[offset++] = edges[j];
    int_qsort(&col_indices[row_ptrs[i]+1], edge_offsets[i+1]-edge_offsets[i]);
  }
  row_ptrs[num_rows] = offset;
  ASSERT(offset == num_nz);

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

// Here's our finite difference implementation of the Jacobian matrix-vector 
// product. 
static void finite_diff_F_Jv(KINSysFn F, void* context, N_Vector u, N_Vector v, N_Vector* work, N_Vector Jv)
{
  static double eps = UNIT_ROUNDOFF;

  // work[1] contains F(u).

  // u + eps*v -> work[2].
  for (int i = 0; i < NV_LOCLENGTH(u); ++i)
    NV_Ith(work[2], i) = NV_Ith(u, i) + eps*NV_Ith(v, i);

  // F(u + eps*v) -> work[3].
  F(work[2], work[3], context);

  // (F(u + eps*v) - F(u)) / eps -> Jv
  for (int i = 0; i < NV_LOCLENGTH(u); ++i)
    NV_Ith(Jv, i) = (NV_Ith(work[3], i) - NV_Ith(work[1], i)) / eps;
}

static void insert_Jv_into_matrix(adj_graph_t* graph, 
                                  adj_graph_coloring_t* coloring, 
                                  int color, 
                                  N_Vector Jv, 
                                  SuperMatrix* J)
{
  int N = NV_LOCLENGTH(Jv);
  NRformat* Jdata = J->Store;
  double* Jij = Jdata->nzval;
  for (int i = 0; i < N; ++i)
  {
    if (NV_Ith(Jv, i) != 0.0)
    {
      int pos = 0, j;
      while (adj_graph_next_edge(graph, i, &pos, &j))
      {
        if (i == j)
        {
          // Diagonal value.
          Jij[Jdata->rowptr[i]] = NV_Ith(Jv, i);
        }
        else if (adj_graph_coloring_has_vertex(coloring, color, j))
        {
          // Off-diagonal value.
          int row_index = Jdata->rowptr[i];
          size_t num_cols = Jdata->colind[Jdata->rowptr[i+1]] - Jdata->colind[row_index];
          int* entry = int_bsearch(&Jdata->colind[row_index+1], num_cols - 1, j);
          ASSERT(entry != NULL);
          Jij[*entry] = NV_Ith(Jv, i);
        }
      }
    }
  }
}

static void compute_F_jacobian(KINSysFn F, 
                               void* context, 
                               N_Vector u, 
                               adj_graph_t* graph, 
                               adj_graph_coloring_t* coloring, 
                               N_Vector* work,
                               SuperMatrix* J)
{
  // We compute the system Jacobian using the method described in 
  // Curtis, Powell, and Reed.
  int N = NV_LOCLENGTH(u);
  N_Vector Jv = N_VClone(u);
  int num_colors = adj_graph_coloring_num_colors(coloring);
  for (int c = 0; c < num_colors; ++c)
  {
    // We construct d, the binary vector corresponding to this color, in work[0].
    memset(work[0], 0, sizeof(double) * N);
    int pos = 0, i;
    while (adj_graph_coloring_next_vertex(coloring, c, &pos, &i))
      NV_Ith(work[0], i) = 1.0;

    // We evaluate F(u) and place it into work[1].
    F(u, work[1], context);

    // Now evaluate the matrix-vector product.
    memset(NV_DATA(Jv), 0, sizeof(double) * N);
    finite_diff_F_Jv(F, context, u, work[0], work, Jv);

    // Copy the components of Jv into their proper locations.
    insert_Jv_into_matrix(graph, coloring, c, Jv, J);
  }
  N_VDestroy(Jv);
}

// Here's our finite difference implementation of the RHS Jacobian 
// matrix-vector product. 
static void finite_diff_rhs_Jv(CVRhsFn rhs, void* context, N_Vector u, double t, N_Vector v, N_Vector* work, N_Vector Jv)
{
  static double eps = UNIT_ROUNDOFF;

  // work[1] contains rhs(u, t).

  // u + eps*v -> work[2].
  for (int i = 0; i < NV_LOCLENGTH(u); ++i)
    NV_Ith(work[2], i) = NV_Ith(u, i) + eps*NV_Ith(v, i);

  // F(u + eps*v, t) -> work[3].
  rhs(t, work[2], work[3], context);

  // (F(u + eps*v) - F(u)) / eps -> Jv
  for (int i = 0; i < NV_LOCLENGTH(u); ++i)
    NV_Ith(Jv, i) = (NV_Ith(work[3], i) - NV_Ith(work[1], i)) / eps;
}

static void compute_rhs_jacobian(CVRhsFn rhs, 
                                 void* context, 
                                 N_Vector u, 
                                 double t, 
                                 adj_graph_t* graph, 
                                 adj_graph_coloring_t* coloring, 
                                 N_Vector* work,
                                 SuperMatrix* J)
{
  // We compute the system Jacobian using the method described in 
  // Curtis, Powell, and Reed.
  int N = NV_LOCLENGTH(u);
  N_Vector Jv = N_VClone(u);
  int num_colors = adj_graph_coloring_num_colors(coloring);
  for (int c = 0; c < num_colors; ++c)
  {
    // We construct d, the binary vector corresponding to this color, in work[0].
    memset(work[0], 0, sizeof(double) * N);
    int pos = 0, i;
    while (adj_graph_coloring_next_vertex(coloring, c, &pos, &i))
      NV_Ith(work[0], i) = 1.0;

    // We evaluate rhs(u, t) and place it in work[1].
    rhs(t, u, work[1], context); 

    // Now evaluate the matrix-vector product.
    memset(NV_DATA(Jv), 0, sizeof(double) * N);
    finite_diff_rhs_Jv(rhs, context, u, t, work[0], work, Jv);

    // Copy the components of Jv into their proper locations.
    insert_Jv_into_matrix(graph, coloring, c, Jv, J);
  }
  N_VDestroy(Jv);
}

void supermatrix_factory_update_jacobian(supermatrix_factory_t* factory, N_Vector u, double t, SuperMatrix* J)
{
  if (factory->F != NULL)
  {
    factory->set_F_time(t, factory->context);
    compute_F_jacobian(factory->F, factory->context, u, factory->graph, factory->coloring, factory->work, J);
  }
  else
  {
    ASSERT(factory->rhs != NULL);
    compute_rhs_jacobian(factory->rhs, factory->context, u, t, factory->graph, factory->coloring, factory->work, J);
  }
}

void supermatrix_free(SuperMatrix* matrix)
{
  Destroy_CompCol_Matrix(matrix);
  free(matrix);
}

