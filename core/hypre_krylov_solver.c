// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <float.h>
#include <dlfcn.h>
#include "core/krylov_solver.h"
#include "core/options.h"
#include "core/timer.h"
#include "core/text_buffer.h"
#include "core/string_utils.h"

//------------------------------------------------------------------------
// This file implements the dynamically-loadable HYPRE Krylov solver.
//------------------------------------------------------------------------

// HYPRE types.
typedef double HYPRE_Real;
typedef double _Complex HYPRE_Complex;
typedef long long HYPRE_Int;
typedef void* HYPRE_Solver;
typedef void* HYPRE_IJMatrix;
typedef void* HYPRE_IJVector;
typedef void* HYPRE_ParCSRMatrix;
typedef void* HYPRE_ParVector;
typedef HYPRE_Int (*HYPRE_PtrToParSolverFcn)(HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector);
#if POLYMEC_HAVE_MPI
typedef MPI_Comm HYPRE_MPI_Comm;
#else
typedef HYPRE_Int HYPRE_MPI_Comm;
#endif

// Here's a table of function pointers for the HYPRE library.
typedef struct
{
  HYPRE_Int (*HYPRE_ParCSRGMRESCreate)(HYPRE_MPI_Comm, HYPRE_Solver*);
  HYPRE_Int (*HYPRE_ParCSRGMRESDestroy)(HYPRE_Solver);
  HYPRE_Int (*HYPRE_ParCSRGMRESSetup)(HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector);
  HYPRE_Int (*HYPRE_ParCSRGMRESSolve)(HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector);
  HYPRE_Int (*HYPRE_ParCSRGMRESSetKDim)(HYPRE_Solver, HYPRE_Int);
  HYPRE_Int (*HYPRE_ParCSRGMRESSetTol)(HYPRE_Solver, HYPRE_Real);
  HYPRE_Int (*HYPRE_ParCSRGMRESSetAbsoluteTol)(HYPRE_Solver, HYPRE_Real);
  HYPRE_Int (*HYPRE_ParCSRGMRESSetMaxIter)(HYPRE_Solver, HYPRE_Int);
  HYPRE_Int (*HYPRE_ParCSRGMRESSetStopCrit)(HYPRE_Solver, HYPRE_Int);
  HYPRE_Int (*HYPRE_ParCSRGMRESSetPrecond)(HYPRE_Solver, HYPRE_PtrToParSolverFcn, HYPRE_PtrToParSolverFcn, HYPRE_Solver);
  HYPRE_Int (*HYPRE_ParCSRGMRESGetPrecond)(HYPRE_Solver, HYPRE_Solver*);
  HYPRE_Int (*HYPRE_ParCSRGMRESGetNumIterations)(HYPRE_Solver, HYPRE_Int*);
  HYPRE_Int (*HYPRE_ParCSRGMRESGetFinalRelativeResidualNorm)(HYPRE_Solver, HYPRE_Real*);

  HYPRE_Int (*HYPRE_ParCSRBiCGSTABCreate)(HYPRE_MPI_Comm, HYPRE_Solver*);
  HYPRE_Int (*HYPRE_ParCSRBiCGSTABDestroy)(HYPRE_Solver);
  HYPRE_Int (*HYPRE_ParCSRBiCGSTABSetup)(HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector);
  HYPRE_Int (*HYPRE_ParCSRBiCGSTABSolve)(HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector);
  HYPRE_Int (*HYPRE_ParCSRBiCGSTABSetTol)(HYPRE_Solver, HYPRE_Real);
  HYPRE_Int (*HYPRE_ParCSRBiCGSTABSetAbsoluteTol)(HYPRE_Solver, HYPRE_Real);
  HYPRE_Int (*HYPRE_ParCSRBiCGSTABSetMaxIter)(HYPRE_Solver, HYPRE_Int);
  HYPRE_Int (*HYPRE_ParCSRBiCGSTABSetStopCrit)(HYPRE_Solver, HYPRE_Int);
  HYPRE_Int (*HYPRE_ParCSRBiCGSTABSetPrecond)(HYPRE_Solver, HYPRE_PtrToParSolverFcn, HYPRE_PtrToParSolverFcn, HYPRE_Solver);
  HYPRE_Int (*HYPRE_ParCSRBiCGSTABGetPrecond)(HYPRE_Solver, HYPRE_Solver*);
  HYPRE_Int (*HYPRE_ParCSRBiCGSTABGetNumIterations)(HYPRE_Solver, HYPRE_Int*);
  HYPRE_Int (*HYPRE_ParCSRBiCGSTABGetFinalRelativeResidualNorm)(HYPRE_Solver, HYPRE_Real*);

  HYPRE_Int (*HYPRE_BoomerAMGCreate)(HYPRE_Solver*);
  HYPRE_Int (*HYPRE_BoomerAMGDestroy)(HYPRE_Solver);
  HYPRE_Int (*HYPRE_BoomerAMGSetup)(HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector);
  HYPRE_Int (*HYPRE_BoomerAMGSolve)(HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector);
  HYPRE_Int (*HYPRE_BoomerAMGGetNumIterations)(HYPRE_Solver, HYPRE_Int*);
  HYPRE_Int (*HYPRE_BoomerAMGGetFinalRelativeResidualNorm)(HYPRE_Solver, HYPRE_Real*);
  HYPRE_Int (*HYPRE_BoomerAMGSetNumFunctions)(HYPRE_Solver, HYPRE_Int);
  HYPRE_Int (*HYPRE_BoomerAMGSetDofFunc)(HYPRE_Solver, HYPRE_Int*);
  HYPRE_Int (*HYPRE_BoomerAMGSetTol)(HYPRE_Solver, HYPRE_Real);
  HYPRE_Int (*HYPRE_BoomerAMGSetMaxIter)(HYPRE_Solver, HYPRE_Int);
  HYPRE_Int (*HYPRE_BoomerAMGSetMaxCoarseSize)(HYPRE_Solver, HYPRE_Int);
  HYPRE_Int (*HYPRE_BoomerAMGSetMinCoarseSize)(HYPRE_Solver, HYPRE_Int);
  HYPRE_Int (*HYPRE_BoomerAMGSetMaxLevels)(HYPRE_Solver, HYPRE_Int);
  HYPRE_Int (*HYPRE_BoomerAMGSetStrongThreshold)(HYPRE_Solver, HYPRE_Real);
  HYPRE_Int (*HYPRE_BoomerAMGSetMaxRowSum)(HYPRE_Solver, HYPRE_Real);
  HYPRE_Int (*HYPRE_BoomerAMGSetCoarsenType)(HYPRE_Solver, HYPRE_Int);
  HYPRE_Int (*HYPRE_BoomerAMGSetMeasureType)(HYPRE_Solver, HYPRE_Int);
  HYPRE_Int (*HYPRE_BoomerAMGSetAggNumLevels)(HYPRE_Solver, HYPRE_Int);
  HYPRE_Int (*HYPRE_BoomerAMGSetNumPaths)(HYPRE_Solver, HYPRE_Int);
  HYPRE_Int (*HYPRE_BoomerAMGSetCGCIts)(HYPRE_Solver, HYPRE_Int);
  HYPRE_Int (*HYPRE_BoomerAMGSetNodal)(HYPRE_Solver, HYPRE_Int);
  HYPRE_Int (*HYPRE_BoomerAMGSetNodalDiag)(HYPRE_Solver, HYPRE_Int);
  HYPRE_Int (*HYPRE_BoomerAMGSetInterpType)(HYPRE_Solver, HYPRE_Int);
  HYPRE_Int (*HYPRE_BoomerAMGSetTruncFactor)(HYPRE_Solver, HYPRE_Real);
  HYPRE_Int (*HYPRE_BoomerAMGSetPMaxElmts)(HYPRE_Solver, HYPRE_Int);

  HYPRE_Int (*HYPRE_EuclidCreate)(HYPRE_MPI_Comm, HYPRE_Solver*);
  HYPRE_Int (*HYPRE_EuclidDestroy)(HYPRE_Solver);
  HYPRE_Int (*HYPRE_EuclidSetup)(HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector);
  HYPRE_Int (*HYPRE_EuclidSetParams)(HYPRE_Solver, int argc, char* argv[]);
  HYPRE_Int (*HYPRE_EuclidSetLevel)(HYPRE_Solver, HYPRE_Int);
  HYPRE_Int (*HYPRE_EuclidSetBJ)(HYPRE_Solver, HYPRE_Int);
  HYPRE_Int (*HYPRE_EuclidSetSparseA)(HYPRE_Solver, HYPRE_Real);
  HYPRE_Int (*HYPRE_EuclidSetRowScale)(HYPRE_Solver, HYPRE_Int);
  HYPRE_Int (*HYPRE_EuclidSetILUT)(HYPRE_Solver, HYPRE_Real);

  HYPRE_Int (*HYPRE_ParCSRParaSailsCreate)(HYPRE_MPI_Comm, HYPRE_Solver*);
  HYPRE_Int (*HYPRE_ParCSRParaSailsDestroy)(HYPRE_Solver); 
  HYPRE_Int (*HYPRE_ParCSRParaSailsSetup)(HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector);
  HYPRE_Int (*HYPRE_ParCSRParaSailsSolve)(HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector);
  HYPRE_Int (*HYPRE_ParCSRParaSailsSetParams)(HYPRE_Solver, HYPRE_Real, HYPRE_Int);
  HYPRE_Int (*HYPRE_ParCSRParaSailsSetFilter)(HYPRE_Solver, HYPRE_Real);
  HYPRE_Int (*HYPRE_ParCSRParaSailsSetSym)(HYPRE_Solver, HYPRE_Int);
  HYPRE_Int (*HYPRE_ParCSRParaSailsSetLoadbal)(HYPRE_Solver, HYPRE_Real);
  HYPRE_Int (*HYPRE_ParCSRParaSailsSetReuse)(HYPRE_Solver, HYPRE_Int);

  HYPRE_Int (*HYPRE_IJMatrixCreate)(HYPRE_MPI_Comm, HYPRE_Int, HYPRE_Int, HYPRE_Int, HYPRE_Int, HYPRE_IJMatrix*);
  HYPRE_Int (*HYPRE_IJMatrixDestroy)(HYPRE_IJMatrix);
  HYPRE_Int (*HYPRE_IJMatrixInitialize)(HYPRE_IJMatrix);
  HYPRE_Int (*HYPRE_IJMatrixSetValues)(HYPRE_IJMatrix, HYPRE_Int, HYPRE_Int*, const HYPRE_Int*, const HYPRE_Int*, const HYPRE_Real*);
  HYPRE_Int (*HYPRE_IJMatrixAddToValues)(HYPRE_IJMatrix, HYPRE_Int, HYPRE_Int*, const HYPRE_Int*, const HYPRE_Int*, const HYPRE_Real*);
  HYPRE_Int (*HYPRE_IJMatrixAssemble)(HYPRE_IJMatrix);
  HYPRE_Int (*HYPRE_IJMatrixGetValues)(HYPRE_IJMatrix, HYPRE_Int, HYPRE_Int*, HYPRE_Int*, HYPRE_Int*, HYPRE_Real*);
  HYPRE_Int (*HYPRE_IJMatrixSetObjectType)(HYPRE_IJMatrix, HYPRE_Int);
  HYPRE_Int (*HYPRE_IJMatrixGetObject)(HYPRE_IJMatrix, void**);
  HYPRE_Int (*HYPRE_IJMatrixSetRowSizes)(HYPRE_IJMatrix, const HYPRE_Int*);
  HYPRE_Int (*HYPRE_IJMatrixSetDiagOffdSizes)(HYPRE_IJMatrix, const HYPRE_Int*, const HYPRE_Int*);

  HYPRE_Int (*HYPRE_IJVectorCreate)(HYPRE_MPI_Comm, HYPRE_Int, HYPRE_Int, HYPRE_IJVector*);
  HYPRE_Int (*HYPRE_IJVectorDestroy)(HYPRE_IJVector);
  HYPRE_Int (*HYPRE_IJVectorInitialize)(HYPRE_IJVector);
  HYPRE_Int (*HYPRE_IJVectorSetValues)(HYPRE_IJVector, HYPRE_Int, const HYPRE_Int*, const HYPRE_Real*);
  HYPRE_Int (*HYPRE_IJVectorAddToValues)(HYPRE_IJVector, HYPRE_Int, const HYPRE_Int*, const HYPRE_Real*);
  HYPRE_Int (*HYPRE_IJVectorAssemble)(HYPRE_IJVector);
  HYPRE_Int (*HYPRE_IJVectorGetValues)(HYPRE_IJVector, HYPRE_Int, HYPRE_Int*, HYPRE_Real*);
  HYPRE_Int (*HYPRE_IJVectorSetObjectType)(HYPRE_IJVector, HYPRE_Int);
  HYPRE_Int (*HYPRE_IJVectorGetObject)(HYPRE_IJVector, void**);
} hypre_methods_table;

typedef struct
{
  void* hypre;
  hypre_methods_table methods;
} hypre_factory_t;

typedef struct
{
  hypre_factory_t* factory;
} hypre_solver_t;

typedef struct
{
  hypre_factory_t* factory;
  HYPRE_IJMatrix* A;
} hypre_matrix_t;

typedef struct
{
  hypre_factory_t* factory;
  HYPRE_IJVector* v;
} hypre_vector_t;

static void hypre_solver_set_tolerances(void* context,
                                        real_t rel_tol,
                                        real_t abs_tol,
                                        real_t div_tol)

{
  hypre_solver_t* solver = context;
//  HYPRE_Real r, a, d;
//  HYPRE_Int iters;
//  solver->factory->methods.KSPGetTolerances(solver->ksp, &r, &a, &d, &iters);
//  solver->factory->methods.KSPSetTolerances(solver->ksp, rel_tol, abs_tol, 
//                                            div_tol, iters);
}

static void hypre_solver_set_max_iterations(void* context,
                                            int max_iters)
{
  hypre_solver_t* solver = context;
//  PetscReal r, a, d;
//  PetscInt iters;
//  solver->factory->methods.KSPGetTolerances(solver->ksp, &r, &a, &d, &iters);
//  solver->factory->methods.KSPSetTolerances(solver->ksp, r, a, d, max_iters);
}

static void hypre_solver_set_operator(void* context,
                                      void* op)
{
  hypre_solver_t* solver = context;
//  Mat A = op;
//  PetscErrorCode err = solver->factory->methods.KSPSetOperators(solver->ksp, A, A);
//  if (err != 0)
//    polymec_error("hypre_solver_set_operator failed!");
}

static bool hypre_solver_solve(void* context,
                               void* x,
                               void* b,
                               real_t* res_norm,
                               int* num_iters)
{
  hypre_solver_t* solver = context;
  return false;
//  Vec X = x;
//  Vec B = b;
//  PetscErrorCode err = solver->factory->methods.KSPSolve(solver->ksp, X, B);
//  if (err != 0)
//    polymec_error("petsc_solver_solve: Call to linear solve failed!");
//  int code;
//  PetscInt niters;
//  solver->factory->methods.KSPGetConvergedReason(solver->ksp, &code);
//  solver->factory->methods.KSPGetIterationNumber(solver->ksp, &niters);
//  *num_iters = (int)niters;
//  solver->factory->methods.KSPGetResidualNorm(solver->ksp, res_norm);
//  return (code > 0);
}

static void hypre_solver_dtor(void* context)
{
  hypre_solver_t* solver = context;
//  solver->factory->methods.KSPDestroy(&solver->ksp);
//  solver->factory = NULL;
  polymec_free(solver);
}

static krylov_solver_t* hypre_factory_solver(void* context,
                                             MPI_Comm comm,
                                             string_string_unordered_map_t* options)
{
  hypre_solver_t* solver = polymec_malloc(sizeof(hypre_solver_t));
  solver->factory = context;
//  solver->factory->methods.KSPCreate(comm, &solver->ksp);

#if 0
  // Default KSP type is GMRES.
  char** type_p = (options != NULL) ? string_string_unordered_map_get(options, "type") : NULL;
  if (type_p == NULL)
    solver->factory->methods.KSPSetType(solver->ksp, "gmres");
  else
    solver->factory->methods.KSPSetType(solver->ksp, *type_p);

  // Handle the preconditioner's options.
  PC pc;
  solver->factory->methods.KSPGetPC(solver->ksp, &pc);
  solver->factory->methods.PCSetFromOptions(pc);

  // Set the thing up.
  solver->factory->methods.KSPSetFromOptions(solver->ksp);
  solver->factory->methods.KSPSetUp(solver->ksp);
#endif

  // Set up the virtual table.
  krylov_solver_vtable vtable = {.set_tolerances = hypre_solver_set_tolerances,
                                 .set_max_iterations = hypre_solver_set_max_iterations,
                                 .set_operator = hypre_solver_set_operator,
                                 .solve = hypre_solver_solve,
                                 .dtor = hypre_solver_dtor};
  return krylov_solver_new("HYPRE", solver, vtable);
}

static void* hypre_matrix_clone(void* context)
{
  hypre_matrix_t* A = context;
  hypre_matrix_t* clone = polymec_malloc(sizeof(hypre_matrix_t));
  clone->factory = A->factory;
//  PetscErrorCode err = A->factory->methods.MatConvert(A->A, MATSAME, MAT_INITIAL_MATRIX, &(clone->A));
//  if (err != 0)
//    polymec_error("hypre_matrix_clone: Error!");
  return clone;
}

static void hypre_matrix_zero(void* context)
{
  hypre_matrix_t* A = context;
//  A->factory->methods.MatZeroEntries(A->A);
}

static void hypre_matrix_scale(void* context, real_t scale_factor)
{
  hypre_matrix_t* A = context;
//  A->factory->methods.MatScale(A->A, scale_factor);
}

static void hypre_matrix_add_identity(void* context, real_t scale_factor)
{
  hypre_matrix_t* A = context;
//  A->factory->methods.MatShift(A->A, scale_factor);
}

static void hypre_matrix_set_diagonal(void* context, void* D)
{
  hypre_matrix_t* A = context;
//  Vec diag = D;
//  A->factory->methods.MatDiagonalSet(A->A, diag, INSERT_VALUES);
}

static void hypre_matrix_add_diagonal(void* context, void* D)
{
  hypre_matrix_t* A = context;
//  Vec diag = D;
//  A->factory->methods.MatDiagonalSet(A->A, diag, INSERT_VALUES);
}

static void hypre_matrix_set_values(void* context, int num_rows,
                                    int* num_columns, int* rows, int* columns,
                                    real_t* values)
{
  hypre_matrix_t* A = context;
  int col_offset = 0;
  for (int r = 0; r < num_rows; ++r)
  {
    int num_cols = num_columns[r];
    HYPRE_Int row = rows[r];
    HYPRE_Int cols[num_cols];
    real_t vals[num_cols];
    for (int c = 0; c < num_cols; ++c, ++col_offset)
    {
      cols[c] = columns[col_offset];
      vals[c] = values[col_offset];
    }
//    A->factory->methods.MatSetValues(A->A, 1, &row, num_cols, cols, vals, insert_mode);
  }
}

static void hypre_matrix_add_values(void* context, int num_rows,
                                    int* num_columns, int* rows, int* columns,
                                    real_t* values)
{
  hypre_matrix_t* A = context;
  int col_offset = 0;
  for (int r = 0; r < num_rows; ++r)
  {
    int num_cols = num_columns[r];
    HYPRE_Int row = rows[r];
    HYPRE_Int cols[num_cols];
    real_t vals[num_cols];
    for (int c = 0; c < num_cols; ++c, ++col_offset)
    {
      cols[c] = columns[col_offset];
      vals[c] = values[col_offset];
    }
//    A->factory->methods.MatSetValues(A->A, 1, &row, num_cols, cols, vals, insert_mode);
  }
}

static void hypre_matrix_start_assembly(void* context)
{
  hypre_matrix_t* A = context;
//  A->factory->methods.MatAssemblyBegin(A->A, MAT_FINAL_ASSEMBLY);
}

static void hypre_matrix_finish_assembly(void* context)
{
  hypre_matrix_t* A = context;
//  A->factory->methods.MatAssemblyEnd(A->A, MAT_FINAL_ASSEMBLY);
}

static void hypre_matrix_get_values(void* context, int num_rows,
                                    int* num_columns, int* rows, int* columns,
                                    real_t* values)
{
  hypre_matrix_t* A = context;
  int col_offset = 0;
  for (int r = 0; r < num_rows; ++r)
  {
    int num_cols = num_columns[r];
    HYPRE_Int row = rows[r];
    HYPRE_Int cols[num_cols];
    for (int c = 0; c < num_cols; ++c, ++col_offset)
      cols[c] = columns[col_offset];
//    A->factory->methods.MatGetValues(A->A, 1, &row, num_cols, cols, &values[col_offset-num_cols]);
  }
}

static void hypre_matrix_dtor(void* context)
{
  hypre_matrix_t* A = context;
//  A->factory->methods.MatDestroy(&(A->A));
  A->factory = NULL;
  polymec_free(A);
}

static krylov_matrix_t* hypre_factory_matrix(void* context,
                                             adj_graph_t* sparsity)
{
  hypre_matrix_t* A = polymec_malloc(sizeof(hypre_matrix_t));
  A->factory = context;
  int N_local = adj_graph_num_vertices(sparsity);
  index_t* vtx_dist = adj_graph_vertex_dist(sparsity);
  MPI_Comm comm = adj_graph_comm(sparsity);
  int nprocs;
  MPI_Comm_size(comm, &nprocs);
  int N_global = (int)vtx_dist[nprocs];
#if 0
  A->factory->methods.MatCreate(comm, &A->A);
  if (nprocs == 1)
    A->factory->methods.MatSetType(A->A, MATSEQAIJ);
  else
    A->factory->methods.MatSetType(A->A, MATMPIAIJ);
  A->factory->methods.MatSetSizes(A->A, N_local, N_local, N_global, N_global);

  // Perform preallocation.
  if (nprocs == 1)
  {
    PetscInt nnz[N_local];
    for (int v = 0; v < N_local; ++v)
      nnz[v] = (PetscInt)(1 + adj_graph_num_edges(sparsity, v));
      A->factory->methods.MatSeqAIJSetPreallocation(A->A, PETSC_DEFAULT, nnz);
  }
  else
  {
    PetscInt d_nnz[N_local], o_nnz[N_local];
    for (int v = 0; v < N_local; ++v)
    {
      d_nnz[v] = 1; // Diagonal entry.
      o_nnz[v] = 0; // Diagonal entry.
      int num_edges = adj_graph_num_edges(sparsity, v);
      int* edges = adj_graph_edges(sparsity, v);
      for (int e = 0; e < num_edges; ++e)
      {
        int edge = edges[e];
        if (edge >= N_local)
          o_nnz[v] += 1;
        else
          d_nnz[v] += 1;
      }
    }
    A->factory->methods.MatMPIAIJSetPreallocation(A->A, PETSC_DEFAULT, d_nnz, PETSC_DEFAULT, o_nnz);
  }
  A->factory->methods.MatSetUp(A->A);
#endif

  // Set up the virtual table.
  krylov_matrix_vtable vtable = {.clone = hypre_matrix_clone,
                                 .zero = hypre_matrix_zero,
                                 .scale = hypre_matrix_scale,
                                 .add_identity = hypre_matrix_add_identity,
                                 .add_diagonal = hypre_matrix_add_diagonal,
                                 .set_diagonal = hypre_matrix_set_diagonal,
                                 .set_values = hypre_matrix_set_values,
                                 .add_values = hypre_matrix_add_values,
                                 .start_assembly = hypre_matrix_start_assembly,
                                 .finish_assembly = hypre_matrix_finish_assembly,
                                 .get_values = hypre_matrix_get_values,
                                 .dtor = hypre_matrix_dtor};
  return krylov_matrix_new(A, vtable, N_local, N_global);
}

static krylov_matrix_t* hypre_factory_block_matrix(void* context,
                                                   adj_graph_t* sparsity,
                                                   int block_size)
{
  hypre_matrix_t* A = polymec_malloc(sizeof(hypre_matrix_t));
  A->factory = context;
  int N_local = adj_graph_num_vertices(sparsity);
  index_t* vtx_dist = adj_graph_vertex_dist(sparsity);
  MPI_Comm comm = adj_graph_comm(sparsity);
  int nprocs;
  MPI_Comm_size(comm, &nprocs);
  int N_global = (int)vtx_dist[nprocs];
#if 0
  A->factory->methods.MatCreate(comm, &A->A);
  if (nprocs == 1)
    A->factory->methods.MatSetType(A->A, MATSEQBAIJ);
  else
    A->factory->methods.MatSetType(A->A, MATMPIBAIJ);
  A->factory->methods.MatSetSizes(A->A, N_local, N_local, N_global, N_global);
  A->factory->methods.MatSetBlockSize(A->A, block_size);

  // Perform preallocation.
  if (nprocs == 1)
  {
    PetscInt nnz[N_local];
    for (int v = 0; v < N_local; ++v)
      nnz[v] = (PetscInt)(1 + adj_graph_num_edges(sparsity, v));
      A->factory->methods.MatSeqBAIJSetPreallocation(A->A, (PetscInt)block_size, PETSC_DEFAULT, nnz);
  }
  else
  {
    PetscInt d_nnz[N_local], o_nnz[N_local];
    for (int v = 0; v < N_local; ++v)
    {
      d_nnz[v] = 1; // Diagonal entry.
      o_nnz[v] = 0; // Diagonal entry.
      int num_edges = adj_graph_num_edges(sparsity, v);
      int* edges = adj_graph_edges(sparsity, v);
      for (int e = 0; e < num_edges; ++e)
      {
        int edge = edges[e];
        if (edge >= N_local)
          o_nnz[v] += 1;
        else
          d_nnz[v] += 1;
      }
    }
    A->factory->methods.MatMPIBAIJSetPreallocation(A->A, (PetscInt)block_size, PETSC_DEFAULT, d_nnz, PETSC_DEFAULT, o_nnz);
  }
  A->factory->methods.MatSetUp(A->A);
#endif

  // Set up the virtual table.
  krylov_matrix_vtable vtable = {.clone = hypre_matrix_clone,
                                 .zero = hypre_matrix_zero,
                                 .scale = hypre_matrix_scale,
                                 .add_identity = hypre_matrix_add_identity,
                                 .add_diagonal = hypre_matrix_add_diagonal,
                                 .set_diagonal = hypre_matrix_set_diagonal,
                                 .set_values = hypre_matrix_set_values,
                                 .add_values = hypre_matrix_add_values,
                                 .start_assembly = hypre_matrix_start_assembly,
                                 .finish_assembly = hypre_matrix_finish_assembly,
                                 .get_values = hypre_matrix_get_values,
                                 .dtor = hypre_matrix_dtor};
  return krylov_matrix_new(A, vtable, N_local, N_global);
}

static void* hypre_vector_clone(void* context)
{
  hypre_vector_t* v = context;
  hypre_vector_t* clone = polymec_malloc(sizeof(hypre_vector_t));
  clone->factory = v->factory;
//  PetscErrorCode err = v->factory->methods.VecDuplicate(v->v, &(clone->v));
//  if (err != 0)
//    polymec_error("hypre_vector_clone: Error in VecDuplicate!");
//  err = v->factory->methods.VecCopy(v->v, clone->v);
//  if (err != 0)
//    polymec_error("hypre_vector_clone: Error in VecCopy!");
  return clone;
}

static void hypre_vector_zero(void* context)
{
  hypre_vector_t* v = context;
//  v->factory->methods.VecZeroEntries(v->v);
}

static void hypre_vector_set_value(void* context, real_t value)
{
  hypre_vector_t* v = context;
//  v->factory->methods.VecSet(v->v, value);
}

static void hypre_vector_scale(void* context, real_t scale_factor)
{
  hypre_vector_t* v = context;
//  v->factory->methods.VecScale(v->v, scale_factor);
}

static void hypre_vector_set_values(void* context, int num_values,
                                    int* indices, real_t* values)
{
  hypre_vector_t* v = context;
//  v->factory->methods.VecSetValues(v->v, num_values, indices, values, INSERT_VALUES);
}

static void hypre_vector_add_values(void* context, int num_values,
                                    int* indices, real_t* values)
{
  hypre_vector_t* v = context;
//  v->factory->methods.VecSetValues(v->v, num_values, indices, values, ADD_VALUES);
}

static void hypre_vector_start_assembly(void* context)
{
  hypre_vector_t* v = context;
//  v->factory->methods.VecAssemblyBegin(v->v);
}

static void hypre_vector_finish_assembly(void* context)
{
  hypre_vector_t* v = context;
//  v->factory->methods.VecAssemblyEnd(v->v);
}

static void hypre_vector_get_values(void* context, int num_values,
                                    int* indices, real_t* values)
{
  hypre_vector_t* v = context;
//  v->factory->methods.VecGetValues(v->v, num_values, indices, values);
}

static real_t hypre_vector_norm(void* context, int p)
{
  hypre_vector_t* v = context;
  real_t norm;
//  NormType norm_type;
//  if (p == 0)
//    norm_type = NORM_INFINITY;
//  else if (p == 1)
//    norm_type = NORM_1;
////  else
//    norm_type = NORM_2;
//  v->factory->methods.VecNorm(v->v, norm_type, &norm);
  return norm;
}

static void hypre_vector_dtor(void* context)
{
  hypre_vector_t* v = context;
//  v->factory->methods.VecDestroy(&(v->v));
  v->factory = NULL;
  polymec_free(v);
}

static krylov_vector_t* hypre_factory_vector(void* context,
                                             MPI_Comm comm,
                                             int N)
{
  hypre_vector_t* v = polymec_malloc(sizeof(hypre_vector_t));
  v->factory = context;
  int nprocs;
  MPI_Comm_size(comm, &nprocs);
//  if (nprocs == 1)
//    v->factory->methods.VecCreateSeq(comm, N, &v->v);
//  else
//    v->factory->methods.VecCreateMPI(comm, N, PETSC_DETERMINE, &v->v);
  HYPRE_Int N_global;
//  v->factory->methods.VecGetSize(v->v, &N_global);
  // Set up the virtual table.
  krylov_vector_vtable vtable = {.clone = hypre_vector_clone,
                                 .zero = hypre_vector_zero,
                                 .set_value = hypre_vector_set_value,
                                 .scale = hypre_vector_scale,
                                 .set_values = hypre_vector_set_values,
                                 .add_values = hypre_vector_add_values,
                                 .start_assembly = hypre_vector_start_assembly,
                                 .finish_assembly = hypre_vector_finish_assembly,
                                 .get_values = hypre_vector_get_values,
                                 .norm = hypre_vector_norm,
                                 .dtor = hypre_vector_dtor};
  return krylov_vector_new(v, vtable, N, N_global);
}

static void hypre_factory_dtor(void* context)
{
  hypre_factory_t* factory = context;
//  if (factory->finalize_hypre)
//  {
//    log_debug("hypre_krylov_factory: Finalizing PETSc.");
//    factory->methods.PetscFinalize();
//  }
  log_debug("hypre_krylov_factory: Closing HYPRE library.");
  dlclose(factory->hypre);
  polymec_free(factory);
}

// Use these to retrieve symbols from dynamically loaded libraries.
#define FETCH_SYMBOL(dylib, symbol_name, function_ptr, fail_label) \
  { \
    void* ptr = dlsym(dylib, symbol_name); \
    if (ptr == NULL) \
    { \
      log_urgent("%s: unable to find %s in dynamic library.", __func__, symbol_name); \
      goto fail_label; \
    } \
    *((void**)&(function_ptr)) = ptr; \
  } 

krylov_factory_t* hypre_krylov_factory(const char* hypre_dir)
{
  hypre_factory_t* factory = polymec_malloc(sizeof(hypre_factory_t));

  // Try to find HYPRE.
  char hypre_path[FILENAME_MAX+1];
  snprintf(hypre_path, FILENAME_MAX, "%s/libHYPRE%s", hypre_dir, SHARED_LIBRARY_SUFFIX);

  // Try to open libPETSc and mine it for symbols.
  log_debug("hypre_krylov_factory: Opening HYPRE library at %s.", hypre_path);
  void* hypre = dlopen(hypre_path, RTLD_NOW);
  if (hypre == NULL)
  {
    char* msg = dlerror();
    polymec_error("hypre_krylov_factory: %s.", msg);
  }

  // If the HYPRE library is parallel, it will contain a reference to MPI symbols, 
  // even if those symbols are not defined within the library itself.
  bool hypre_is_parallel = (dlsym(hypre, "MPI_Abort") != NULL);
#if POLYMEC_HAVE_MPI
  if (!hypre_is_parallel)
  {
    log_urgent("hypre_krylov_factory: Polymec is configured to use MPI, but HYPRE is not.\n"
               "  HYPRE must be built using -DHYPRE_SEQUENTIAL=OFF."); 
    goto failure;
  }
#else
  if (hypre_is_parallel)
  {
    log_urgent("hypre_krylov_factory: Polymec is configured serially, but HYPRE is not.\n"
               "  HYPRE must be built using -DHYPRE_SEQUENTIAL=ON."); 
    goto failure;
  }
#endif
  log_debug("hypre_krylov_factory: Succeeded.");

  // HYPRE defines a function called hypre_printf if it is configured for "big" (64-bit) integers.
  bool hypre_uses_64bit_indices = (dlsym(hypre, "hypre_printf") != NULL);
  if (sizeof(index_t) == sizeof(int64_t))
  {
    if (!hypre_uses_64bit_indices)
    {
      log_urgent("hypre_krylov_factory: Since polymec is configured for 64-bit indices,\n"
                 "  HYPRE must be built using -DHYPRE_BIGINT=ON.");
      goto failure;
    }
  }
  else
  {
    if (hypre_uses_64bit_indices)
    {
      log_urgent("hypre_krylov_factory: Since polymec is configured for 32-bit indices,\n"
                 "  HYPRE must be built using -DHYPRE_BIGINT=OFF.");
      goto failure;
    }
  }

  // Get the symbols.
#define FETCH_HYPRE_SYMBOL(symbol_name) \
  FETCH_SYMBOL(hypre, #symbol_name, factory->methods.symbol_name, failure);

  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRGMRESCreate);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRGMRESDestroy);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRGMRESSetup);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRGMRESSolve);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRGMRESSetKDim);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRGMRESSetTol);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRGMRESSetAbsoluteTol);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRGMRESSetMaxIter);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRGMRESSetStopCrit);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRGMRESSetPrecond);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRGMRESGetPrecond);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRGMRESGetNumIterations);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRGMRESGetFinalRelativeResidualNorm);

  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRBiCGSTABCreate);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRBiCGSTABDestroy);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRBiCGSTABSetup);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRBiCGSTABSolve);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRBiCGSTABSetTol);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRBiCGSTABSetAbsoluteTol);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRBiCGSTABSetMaxIter);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRBiCGSTABSetStopCrit);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRBiCGSTABSetPrecond);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRBiCGSTABGetPrecond);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRBiCGSTABGetNumIterations);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRBiCGSTABGetFinalRelativeResidualNorm);

  FETCH_HYPRE_SYMBOL(HYPRE_BoomerAMGCreate);
  FETCH_HYPRE_SYMBOL(HYPRE_BoomerAMGDestroy);
  FETCH_HYPRE_SYMBOL(HYPRE_BoomerAMGSetup);
  FETCH_HYPRE_SYMBOL(HYPRE_BoomerAMGSolve);
  FETCH_HYPRE_SYMBOL(HYPRE_BoomerAMGGetNumIterations);
  FETCH_HYPRE_SYMBOL(HYPRE_BoomerAMGGetFinalRelativeResidualNorm);
  FETCH_HYPRE_SYMBOL(HYPRE_BoomerAMGSetNumFunctions);
  FETCH_HYPRE_SYMBOL(HYPRE_BoomerAMGSetDofFunc);
  FETCH_HYPRE_SYMBOL(HYPRE_BoomerAMGSetTol);
  FETCH_HYPRE_SYMBOL(HYPRE_BoomerAMGSetMaxIter);
  FETCH_HYPRE_SYMBOL(HYPRE_BoomerAMGSetMaxCoarseSize);
  FETCH_HYPRE_SYMBOL(HYPRE_BoomerAMGSetMinCoarseSize);
  FETCH_HYPRE_SYMBOL(HYPRE_BoomerAMGSetMaxLevels);
  FETCH_HYPRE_SYMBOL(HYPRE_BoomerAMGSetStrongThreshold);
  FETCH_HYPRE_SYMBOL(HYPRE_BoomerAMGSetMaxRowSum);
  FETCH_HYPRE_SYMBOL(HYPRE_BoomerAMGSetCoarsenType);
  FETCH_HYPRE_SYMBOL(HYPRE_BoomerAMGSetMeasureType);
  FETCH_HYPRE_SYMBOL(HYPRE_BoomerAMGSetAggNumLevels);
  FETCH_HYPRE_SYMBOL(HYPRE_BoomerAMGSetNumPaths);
  FETCH_HYPRE_SYMBOL(HYPRE_BoomerAMGSetCGCIts);
  FETCH_HYPRE_SYMBOL(HYPRE_BoomerAMGSetNodal);
  FETCH_HYPRE_SYMBOL(HYPRE_BoomerAMGSetNodalDiag);
  FETCH_HYPRE_SYMBOL(HYPRE_BoomerAMGSetInterpType);
  FETCH_HYPRE_SYMBOL(HYPRE_BoomerAMGSetTruncFactor);
  FETCH_HYPRE_SYMBOL(HYPRE_BoomerAMGSetPMaxElmts);

  FETCH_HYPRE_SYMBOL(HYPRE_EuclidCreate);
  FETCH_HYPRE_SYMBOL(HYPRE_EuclidDestroy);
  FETCH_HYPRE_SYMBOL(HYPRE_EuclidSetup);
  FETCH_HYPRE_SYMBOL(HYPRE_EuclidSetParams);
  FETCH_HYPRE_SYMBOL(HYPRE_EuclidSetLevel);
  FETCH_HYPRE_SYMBOL(HYPRE_EuclidSetBJ);
  FETCH_HYPRE_SYMBOL(HYPRE_EuclidSetSparseA);
  FETCH_HYPRE_SYMBOL(HYPRE_EuclidSetRowScale);
  FETCH_HYPRE_SYMBOL(HYPRE_EuclidSetILUT);

  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRParaSailsCreate);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRParaSailsDestroy);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRParaSailsSetup);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRParaSailsSolve);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRParaSailsSetParams);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRParaSailsSetFilter);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRParaSailsSetSym);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRParaSailsSetLoadbal);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRParaSailsSetReuse);

  FETCH_HYPRE_SYMBOL(HYPRE_IJMatrixCreate);
  FETCH_HYPRE_SYMBOL(HYPRE_IJMatrixDestroy);
  FETCH_HYPRE_SYMBOL(HYPRE_IJMatrixInitialize);
  FETCH_HYPRE_SYMBOL(HYPRE_IJMatrixSetValues);
  FETCH_HYPRE_SYMBOL(HYPRE_IJMatrixAddToValues);
  FETCH_HYPRE_SYMBOL(HYPRE_IJMatrixAssemble);
  FETCH_HYPRE_SYMBOL(HYPRE_IJMatrixGetValues);
  FETCH_HYPRE_SYMBOL(HYPRE_IJMatrixSetObjectType);
  FETCH_HYPRE_SYMBOL(HYPRE_IJMatrixGetObject);
  FETCH_HYPRE_SYMBOL(HYPRE_IJMatrixSetRowSizes);
  FETCH_HYPRE_SYMBOL(HYPRE_IJMatrixSetDiagOffdSizes);

  FETCH_HYPRE_SYMBOL(HYPRE_IJVectorCreate);
  FETCH_HYPRE_SYMBOL(HYPRE_IJVectorDestroy);
  FETCH_HYPRE_SYMBOL(HYPRE_IJVectorInitialize);
  FETCH_HYPRE_SYMBOL(HYPRE_IJVectorSetValues);
  FETCH_HYPRE_SYMBOL(HYPRE_IJVectorAddToValues);
  FETCH_HYPRE_SYMBOL(HYPRE_IJVectorAssemble);
  FETCH_HYPRE_SYMBOL(HYPRE_IJVectorGetValues);
  FETCH_HYPRE_SYMBOL(HYPRE_IJVectorSetObjectType);
  FETCH_HYPRE_SYMBOL(HYPRE_IJVectorGetObject);
#undef FETCH_HYPRE_SYMBOL

  log_debug("hypre_krylov_factory: Got HYPRE symbols.");

  // Stash the library.
  factory->hypre = hypre; 

  // Construct the factory.
  krylov_factory_vtable vtable = {.solver = hypre_factory_solver,
                                  .matrix = hypre_factory_matrix,
                                  .block_matrix = hypre_factory_block_matrix,
                                  .vector = hypre_factory_vector,
                                  .dtor = hypre_factory_dtor};
  return krylov_factory_new("2.10", factory, vtable);

failure:
  dlclose(hypre);
  polymec_free(factory);
  return NULL;
}

