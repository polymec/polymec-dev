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
#include "kinsol/kinsol.h"
#include "kinsol/kinsol_spgmr.h"
#include "kinsol/kinsol_spbcgs.h"
#include "kinsol/kinsol_sptfqmr.h"

typedef struct 
{
  void* context;
  int block_size;
  nonlinear_solver_dtor dtor;
  mf_nonlinear_solver_type_t type;
  adj_graph_t* graph;
  void* kinsol;
} mf_nonlinear_solver_t;

static void compute_precond_block_jacobian(void* context, 
                                           nonlinear_solver_eval_func eval,
                                           double t, 
                                           int block_row, 
                                           int* block_cols, 
                                           int num_block_cols, 
                                           int N,
                                           double* X, 
                                           double** jacobian_cols)
{
  mf_nonlinear_solver_t* mf_solver = context;
  static const double epsilon = 1e-8;

  int block_size = mf_solver->block_size;
  for (int row = block_size*block_row; row < block_size*(block_row + 1); ++row)
  {
    double F0 = eval(mf_solver->context, t, row, N, X);
    // Diagonal terms.
    for (int col = block_size*row; col < block_size*(row+1); ++col)
    {
      double Xj = X[col];
      double dXj = (Xj == 0.0) ? epsilon : epsilon * Xj;
      X[col] = Xj + dXj;
      double F1 = eval(mf_solver->context, t, row, N, X);
      X[col] = Xj;
      jacobian_cols[row][col] = (F1 - F0) / dXj;
    }
    // Off-diagonal terms (local domain only).
    for (int j = 0; j < num_block_cols; ++j)
    {
      for (int col = block_size*block_cols[j]; col < block_size*(block_cols[j]+1); ++col)
      {
        double Xj = X[col];
        double dXj = (Xj == 0.0) ? epsilon : epsilon * Xj;
        X[col] = Xj + dXj;
        double F1 = eval(mf_solver->context, t, row, N, X);
        X[col] = Xj;
        jacobian_cols[row][j+1] = (F1 - F0) / dXj;
      }
    }
  }
}

static void mf_nonlinear_solver_solve(void* context, double t, int N, double* X)
{
  mf_nonlinear_solver_t* mf_solver = context;
}

static void mf_nonlinear_solver_free(void* context)
{
  mf_nonlinear_solver_t* mf_solver = context;
  if ((mf_solver->dtor != NULL) && (mf_solver->context != NULL))
    mf_solver->dtor(mf_solver->context);
  KINFree(&mf_solver->kinsol);
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
  ASSERT(num_equations_per_site > 0);
  mf_solver->context = context;
  mf_solver->block_size = num_equations_per_site;
  mf_solver->graph = graph; // borrowed.
  mf_solver->dtor = dtor;
  mf_solver->type = type;
  nonlinear_solver_vtable vtable = {.eval = eval, 
                                    .solve = mf_nonlinear_solver_solve,
                                    .dtor = mf_nonlinear_solver_free};
  int num_sites = adj_graph_num_vertices(graph);
  int N = num_sites * num_equations_per_site;

  // Set up the kinsol solver object.
  mf_solver->kinsol = KINCreate();
  // FIXME

  return nonlinear_solver_new(name, mf_solver, vtable, N);
}

