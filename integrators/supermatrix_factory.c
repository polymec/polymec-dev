// Copyright (c) 2012-2014, Jeffrey N. Johnson
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

#include "core/array_utils.h"
#include "core/sundials_helpers.h"
#include "integrators/supermatrix_factory.h"
#include "slu_util.h"

struct supermatrix_factory_t 
{
  adj_graph_t* graph;
  adj_graph_coloring_t* coloring;
  KINSysFn F;
  void (*set_F_time)(void* context, real_t t);
  CVRhsFn rhs;
  void* context;

  // Work vectors.
  N_Vector* work;
};

supermatrix_factory_t* supermatrix_factory_from_sys_func(adj_graph_t* graph,
                                                         KINSysFn F,
                                                         void (*set_F_time)(void*, real_t),
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
  real_t* mat_zeros = SUPERLU_MALLOC(sizeof(real_t) * num_nz);
  memset(mat_zeros, 0, sizeof(real_t) * num_nz);

  // Hand over these resources to create the Supermatrix.
  dCreate_CompRow_Matrix(A, num_rows, num_rows, num_nz, 
                         mat_zeros, col_indices, row_ptrs, 
                         SLU_NR, SLU_D, SLU_GE);

  return A;
}

SuperMatrix* supermatrix_factory_vector(supermatrix_factory_t* factory,
                                        int num_rhs)
{
  SuperMatrix* B = malloc(sizeof(SuperMatrix));
  
  // Fetch numrows information from the graph.
  int num_rows = adj_graph_num_vertices(factory->graph);
  real_t* b_mem;
  if ( !(b_mem = SUPERLU_MALLOC(sizeof(real_t) * num_rows * num_rhs)) ) {
    ABORT("Malloc fails for SuperLU b/rhs vector.");
  }
  dCreate_Dense_Matrix(B, num_rows, num_rhs,
                       b_mem, num_rows,
                       SLU_DN, SLU_D, SLU_GE);
  return B;
}

SuperMatrix* supermatrix_factory_jacobian(supermatrix_factory_t* factory, N_Vector x, real_t t)
{
  SuperMatrix* J = supermatrix_factory_matrix(factory);
  supermatrix_factory_update_jacobian(factory, x, t, J);
  return J;
}

// Here's our finite difference implementation of the Jacobian matrix-vector 
// product. 
static void finite_diff_F_Jv(KINSysFn F, void* context, N_Vector x, N_Vector v, N_Vector* work, N_Vector Jv)
{
  real_t eps = sqrt(UNIT_ROUNDOFF);

  // work[0] == v
  // work[1] contains F(x).
  // work[2] == u + eps*v
  // work[3] == F(x + eps*v)

  // u + eps*v -> work[2].
  for (int i = 0; i < NV_LOCLENGTH(x); ++i)
    NV_Ith(work[2], i) = NV_Ith(x, i) + eps*NV_Ith(v, i);

  // F(x + eps*v) -> work[3].
  F(work[2], work[3], context);

  // (F(x + eps*v) - F(x)) / eps -> Jv
  for (int i = 0; i < NV_LOCLENGTH(x); ++i)
    NV_Ith(Jv, i) = (NV_Ith(work[3], i) - NV_Ith(work[1], i)) / eps;
}

static void insert_Jv_into_matrix(adj_graph_t* graph, 
                                  adj_graph_coloring_t* coloring, 
                                  int color, 
                                  N_Vector Jv, 
                                  SuperMatrix* J)
{
  NRformat* Jdata = J->Store;
  real_t* Jij = Jdata->nzval;
  int pos = 0, i;
  while (adj_graph_coloring_next_vertex(coloring, color, &pos, &i))
  {
    if (NV_Ith(Jv, i) != 0.0)
    {
      // Fill in the diagonal element.
      Jij[Jdata->rowptr[i]] = NV_Ith(Jv, i);

      int pos = 0, j;
      while (adj_graph_next_edge(graph, i, &pos, &j))
      {
        // Off-diagonal value.
        int row_index = Jdata->rowptr[i];
        size_t num_cols = Jdata->rowptr[i+1] - row_index;
        int* entry = int_bsearch(&Jdata->colind[row_index+1], num_cols - 1, j);
        ASSERT(entry != NULL);
        size_t offset = entry - &Jdata->colind[row_index];
        Jij[Jdata->rowptr[i] + offset] = NV_Ith(Jv, j);
      }
    }
  }
}

static void compute_F_jacobian(KINSysFn F, 
                               void* context, 
                               N_Vector x, 
                               adj_graph_t* graph, 
                               adj_graph_coloring_t* coloring, 
                               N_Vector* work,
                               SuperMatrix* J)
{
  // We compute the system Jacobian using the method described in 
  // Curtis, Powell, and Reed.
  int N = NV_LOCLENGTH(x);
  N_Vector Jv = N_VClone(x);
  int num_colors = adj_graph_coloring_num_colors(coloring);
  for (int c = 0; c < num_colors; ++c)
  {
    // We construct d, the binary vector corresponding to this color, in work[0].
    memset(NV_DATA(work[0]), 0, sizeof(real_t) * N);
    int pos = 0, i;
    while (adj_graph_coloring_next_vertex(coloring, c, &pos, &i))
    {
      NV_Ith(work[0], i) = 1.0;
    }

    // We evaluate F(x) and place it into work[1].
    F(x, work[1], context);

    // Now evaluate the matrix-vector product.
    memset(NV_DATA(Jv), 0, sizeof(real_t) * N);
    finite_diff_F_Jv(F, context, x, work[0], work, Jv);

    // Copy the components of Jv into their proper locations.
    insert_Jv_into_matrix(graph, coloring, c, Jv, J);
  }
  N_VDestroy(Jv);
}

// Here's our finite difference implementation of the RHS Jacobian 
// matrix-vector product. 
static void finite_diff_rhs_Jv(CVRhsFn rhs, void* context, N_Vector x, real_t t, N_Vector v, N_Vector* work, N_Vector Jv)
{
  real_t eps = sqrt(UNIT_ROUNDOFF);

  // work[1] contains rhs(x, t).

  // u + eps*v -> work[2].
  for (int i = 0; i < NV_LOCLENGTH(x); ++i)
    NV_Ith(work[2], i) = NV_Ith(x, i) + eps*NV_Ith(v, i);

  // F(x + eps*v, t) -> work[3].
  rhs(t, work[2], work[3], context);

  // (F(x + eps*v) - F(x)) / eps -> Jv
  for (int i = 0; i < NV_LOCLENGTH(x); ++i)
    NV_Ith(Jv, i) = (NV_Ith(work[3], i) - NV_Ith(work[1], i)) / eps;
}

static void compute_rhs_jacobian(CVRhsFn rhs, 
                                 void* context, 
                                 N_Vector x, 
                                 real_t t, 
                                 adj_graph_t* graph, 
                                 adj_graph_coloring_t* coloring, 
                                 N_Vector* work,
                                 SuperMatrix* J)
{
  // We compute the system Jacobian using the method described in 
  // Curtis, Powell, and Reed.
  int N = NV_LOCLENGTH(x);
  N_Vector Jv = N_VClone(x);
  int num_colors = adj_graph_coloring_num_colors(coloring);
  for (int c = 0; c < num_colors; ++c)
  {
    // We construct d, the binary vector corresponding to this color, in work[0].
    memset(NV_DATA(work[0]), 0, sizeof(real_t) * N);
    int pos = 0, i;
    while (adj_graph_coloring_next_vertex(coloring, c, &pos, &i))
      NV_Ith(work[0], i) = 1.0;

    // We evaluate rhs(x, t) and place it in work[1].
    rhs(t, x, work[1], context); 

    // Now evaluate the matrix-vector product.
    memset(NV_DATA(Jv), 0, sizeof(real_t) * N);
    finite_diff_rhs_Jv(rhs, context, x, t, work[0], work, Jv);

    // Copy the components of Jv into their proper locations.
    insert_Jv_into_matrix(graph, coloring, c, Jv, J);
  }
  N_VDestroy(Jv);
}

void supermatrix_factory_update_jacobian(supermatrix_factory_t* factory, N_Vector x, real_t t, SuperMatrix* J)
{
  if (factory->F != NULL)
  {
    if (factory->set_F_time != NULL) {
      factory->set_F_time(factory->context, t);
    }
    compute_F_jacobian(factory->F, factory->context, x, factory->graph, factory->coloring, factory->work, J);
  }
  else
  {
    ASSERT(factory->rhs != NULL);
    compute_rhs_jacobian(factory->rhs, factory->context, x, t, factory->graph, factory->coloring, factory->work, J);
  }
}

void supermatrix_free(SuperMatrix* matrix)
{
  switch(matrix->Stype) {
  case SLU_DN:
    Destroy_Dense_Matrix(matrix);
    break;
  case SLU_NR:
    Destroy_CompRow_Matrix(matrix);
    break;
  case SLU_NC:
    Destroy_CompCol_Matrix(matrix);
    break;
  default:
    printf("ERROR: unknown matrix type passed to supermatrix_free!");
    break;
  }
  free(matrix);
}

