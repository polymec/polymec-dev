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

#include "integrators/mf_nonlinear_solver.h"

// We use KINSOL for doing the matrix-free nonlinear solve.
#include "kinsol/kinsol.h"
#include "kinsol/kinsol_spgmr.h"
#include "kinsol/kinsol_spbcgs.h"
#include "kinsol/kinsol_sptfqmr.h"

// We use a serial version of SuperLU to do preconditioning.
#include "slu_ddefs.h"
#include "supermatrix.h"
#include "slu_util.h"

typedef struct 
{
  // Parallel stuff.
  int rank, nprocs;

  // Solver stuff.
  void* context;
  int block_size;
  nonlinear_solver_dtor dtor;
  mf_nonlinear_solver_type_t solver_type;
  adj_graph_t* graph;
  void* kinsol;

  // Preconditioning stuff.
  SuperMatrix precond_mat, precond_rhs, precond_L, precond_U;
  int *precond_rperm, *precond_cperm;
} mf_nonlinear_solver_t;

static void set_up_preconditioner(mf_nonlinear_solver_t* mf_solver)
{
  adj_graph_t* graph = mf_solver->graph;
  int* edge_offsets = adj_graph_edge_offsets(graph);
  int block_size = mf_solver->block_size;
  int num_rows = adj_graph_num_vertices(graph);
  int num_nz = block_size * block_size * edge_offsets[mf_solver->rank+1];

  // Translate the sparsity information in the graph to row indices and 
  // column pointers.
  int* row_indices = intMalloc(num_nz);
  int* col_ptrs = intMalloc(num_rows + 1);

  // Make some zeros to feed to the preconditioner system for initialization.
  double* mat_zeros = doubleMalloc(num_nz);
  memset(mat_zeros, 0, sizeof(double) * num_nz);
  double* rhs_zeros = doubleMalloc(num_rows);
  memset(rhs_zeros, 0, sizeof(double) * num_rows);

  // Create the preconditioner matrix.
  dCreate_CompCol_Matrix(&mf_solver->precond_mat, num_rows, num_rows, num_nz, 
                         mat_zeros, row_indices, col_ptrs, SLU_NC, SLU_D, SLU_GE);

  // Create the preconditioner right-hand side vector (an Nx1 dense matrix).
  dCreate_Dense_Matrix(&mf_solver->precond_rhs, num_rows, 1, 
                       rhs_zeros, num_rows, SLU_DN, SLU_D, SLU_GE);

  // Create permutation vectors for rows and columns for the preconditioner 
  // system.
  mf_solver->precond_rperm = intMalloc(num_rows);
  mf_solver->precond_cperm = intMalloc(num_rows);
}

static void free_preconditioner(mf_nonlinear_solver_t* mf_solver)
{
  SUPERLU_FREE(mf_solver->precond_cperm);
  SUPERLU_FREE(mf_solver->precond_rperm);
//  Destroy_CompCol_Matrix(&mf_solver->precond_U);
//  Destroy_SuperNode_Matrix(&mf_solver->precond_L);
  Destroy_SuperMatrix_Store(&mf_solver->precond_rhs);
  Destroy_CompCol_Matrix(&mf_solver->precond_mat);
}

// This function computes a block row in the Jacobian preconditioner matrix 
// and inserts it into the given sparse matrix.
static void compute_precond_jacobian_block_row(void* context, 
                                               nonlinear_solver_eval_func eval,
                                               double t, 
                                               int block_row, 
                                               int* block_cols, 
                                               int num_block_cols, 
                                               int N,
                                               double* X, 
                                               SuperMatrix* jacobian)
{
  mf_nonlinear_solver_t* mf_solver = context;
  static const double epsilon = 1e-8;

  // We will temporarily store the values of the Jacobian here.
  int block_size = mf_solver->block_size;
  double jacobian_cols[block_size][block_size*(1+num_block_cols)];

  for (int i = 0; i < block_size; ++i)
  {
    int row = block_size*block_row + i;
    double F0 = eval(mf_solver->context, t, row, N, X);
    // Diagonal terms.
    for (int j = 0; j < block_size; ++j)
    {
      int col = block_size*row + j;
      double Xj = X[col];
      double dXj = (Xj == 0.0) ? epsilon : epsilon * Xj;
      X[col] = Xj + dXj;
      double F1 = eval(mf_solver->context, t, row, N, X);
      X[col] = Xj;
      jacobian_cols[i][j] = (F1 - F0) / dXj;
    }
    // Off-diagonal terms (local domain only).
    for (int j = 0; j < num_block_cols; ++j)
    {
      for (int k = 0; k < block_size; ++k)
      {
        int col = block_size*block_cols[j] + k;
        double Xj = X[col];
        double dXj = (Xj == 0.0) ? epsilon : epsilon * Xj;
        X[col] = Xj + dXj;
        double F1 = eval(mf_solver->context, t, row, N, X);
        X[col] = Xj;
        jacobian_cols[i][block_size*(j+1)+k] = (F1 - F0) / dXj;
      }
    }
  }

  // Grab the column-compressed matrix storage.
  struct NCformat* cc_mat = jacobian->Store;

  // Inject the values into the Jacobian matrix.
}

static void mf_nonlinear_solver_solve(void* context, double t, int N, double* X)
{
  mf_nonlinear_solver_t* mf_solver = context;
}

static void mf_nonlinear_solver_free(void* context)
{
  mf_nonlinear_solver_t* mf_solver = context;
  free_preconditioner(mf_solver);
  KINFree(&mf_solver->kinsol);
  if ((mf_solver->dtor != NULL) && (mf_solver->context != NULL))
    mf_solver->dtor(mf_solver->context);
  free(mf_solver);
}

nonlinear_solver_t* mf_nonlinear_solver_new(const char* name,
                                            void* context,
                                            nonlinear_solver_eval_func eval,
                                            nonlinear_solver_dtor dtor,
                                            adj_graph_t* graph,
                                            int num_equations_per_site,
                                            mf_nonlinear_solver_type_t type)
{
  ASSERT(eval != NULL);
  ASSERT(graph != NULL);
  ASSERT(num_equations_per_site > 0);

  mf_nonlinear_solver_t* mf_solver = malloc(sizeof(mf_nonlinear_solver_t));

  MPI_Comm comm = adj_graph_comm(graph);
  MPI_Comm_size(comm, &mf_solver->nprocs);
  MPI_Comm_rank(comm, &mf_solver->rank);

  mf_solver->context = context;
  mf_solver->block_size = num_equations_per_site;
  mf_solver->graph = graph; // borrowed.
  mf_solver->dtor = dtor;
  mf_solver->solver_type = type;

  nonlinear_solver_vtable vtable = {.eval = eval, 
                                    .solve = mf_nonlinear_solver_solve,
                                    .dtor = mf_nonlinear_solver_free};
  int num_sites = adj_graph_num_vertices(graph);
  int N = num_sites * num_equations_per_site;

  // Set up the kinsol solver object.
  mf_solver->kinsol = KINCreate();
  // FIXME

  // Set up the preconditioner matrix and right hand side vector.
  // FIXME: We currently assume that the sparsity graph is static.
  set_up_preconditioner(mf_solver);

  return nonlinear_solver_new(name, mf_solver, vtable, N);
}

