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
// This file implements the dynamically-loadable PETSc Krylov solver.
//------------------------------------------------------------------------

// Here's a table of function pointers for the PETSc library.
typedef real_t PetscScalar;
typedef real_t PetscReal;
typedef int PetscMPIInt;
typedef index_t PetscInt; 
typedef enum { PETSC_FALSE,PETSC_TRUE } PetscBool;
typedef int PetscErrorCode;
typedef void* KSP;
typedef const char* KSPType;
typedef void* PC;
typedef const char* PCType;
typedef void* Mat;
typedef const char* MatType;
typedef void* Vec;
typedef enum {NORM_1=0,NORM_2=1,NORM_FROBENIUS=2,NORM_INFINITY=3,NORM_1_AND_2=4} NormType;
typedef enum {MAT_FLUSH_ASSEMBLY=1,MAT_FINAL_ASSEMBLY=0} MatAssemblyType;
typedef enum {INSERT_VALUES=1,ADD_VALUES=0} InsertMode;
typedef enum {MAT_INITIAL_MATRIX,MAT_REUSE_MATRIX,MAT_IGNORE_MATRIX} MatReuse;
const char* MATSEQAIJ = "seqaij";
const char* MATMPIAIJ = "mpiaij";
const char* MATSEQBAIJ = "seqbaij";
const char* MATMPIBAIJ = "mpibaij";
const char* MATSAME = "same";
#define PETSC_DEFAULT -2
#define PETSC_DECIDE -1
#define PETSC_DETERMINE PETSC_DECIDE
typedef struct
{
  PetscErrorCode (*PetscInitialize)(int*, char***, const char[], const char[]);
  PetscErrorCode (*PetscInitialized)(PetscBool*);
  PetscErrorCode (*PetscFinalize)();

  PetscErrorCode (*KSPCreate)(MPI_Comm,KSP *);
  PetscErrorCode (*KSPSetType)(KSP,KSPType);
  PetscErrorCode (*KSPSetUp)(KSP);
  PetscErrorCode (*KSPGetPC)(KSP,PC*);
  PetscErrorCode (*KSPSetFromOptions)(KSP);
  PetscErrorCode (*KSPSetTolerances)(KSP,PetscReal,PetscReal,PetscReal,PetscInt);
  PetscErrorCode (*KSPGetTolerances)(KSP,PetscReal*,PetscReal*,PetscReal*,PetscInt*);
  PetscErrorCode (*KSPSetOperators)(KSP,Mat,Mat);
  PetscErrorCode (*KSPSolve)(KSP,Vec,Vec);
  PetscErrorCode (*KSPGetConvergedReason)(KSP,int*);
  PetscErrorCode (*KSPGetIterationNumber)(KSP,PetscInt*);
  PetscErrorCode (*KSPGetResidualNorm)(KSP,PetscReal*);
  PetscErrorCode (*KSPDestroy)(KSP*);

  PetscErrorCode (*PCCreate)(MPI_Comm,PC*);
  PetscErrorCode (*PCSetType)(PC, const char*);
  PetscErrorCode (*PCSetFromOptions)(PC);

  PetscErrorCode (*MatCreate)(MPI_Comm,Mat*);
  PetscErrorCode (*MatConvert)(Mat, MatType, MatReuse, Mat*);
  PetscErrorCode (*MatSetType)(Mat, MatType);
  PetscErrorCode (*MatSetSizes)(Mat, PetscInt m, PetscInt n, PetscInt M, PetscInt N);
  PetscErrorCode (*MatSeqAIJSetPreallocation)(Mat, PetscInt, PetscInt[]);
  PetscErrorCode (*MatMPIAIJSetPreallocation)(Mat, PetscInt, PetscInt[], PetscInt, PetscInt[]);
  PetscErrorCode (*MatSeqBAIJSetPreallocation)(Mat, PetscInt, PetscInt, PetscInt[]);
  PetscErrorCode (*MatMPIBAIJSetPreallocation)(Mat, PetscInt, PetscInt, PetscInt[], PetscInt, PetscInt[]);
  PetscErrorCode (*MatSetBlockSize)(Mat, PetscInt);
  PetscErrorCode (*MatSetUp)(Mat);
  PetscErrorCode (*MatDestroy)(Mat*);
  PetscErrorCode (*MatScale)(Mat,PetscScalar);
  PetscErrorCode (*MatShift)(Mat,PetscScalar);
  PetscErrorCode (*MatDiagonalSet)(Mat,Vec,InsertMode);
  PetscErrorCode (*MatZeroEntries)(Mat);
  PetscErrorCode (*MatGetSize)(Mat,PetscInt*,PetscInt*);
  PetscErrorCode (*MatGetLocalSize)(Mat,PetscInt*,PetscInt*);
  PetscErrorCode (*MatSetValues)(Mat,PetscInt,const PetscInt[],PetscInt,const PetscInt[],const PetscScalar[],InsertMode);
  PetscErrorCode (*MatAssemblyBegin)(Mat,MatAssemblyType);
  PetscErrorCode (*MatAssemblyEnd)(Mat,MatAssemblyType);
  PetscErrorCode (*MatGetValues)(Mat,PetscInt,const PetscInt[],PetscInt,const PetscInt[],PetscScalar[]);

  PetscErrorCode (*VecCreateSeq)(MPI_Comm,PetscInt,Vec*);
  PetscErrorCode (*VecCreateMPI)(MPI_Comm,PetscInt,PetscInt,Vec*);
  PetscErrorCode (*VecDuplicate)(Vec,Vec*);
  PetscErrorCode (*VecCopy)(Vec,Vec);
  PetscErrorCode (*VecSetUp)(Vec);
  PetscErrorCode (*VecDestroy)(Vec*);
  PetscErrorCode (*VecGetSize)(Vec, PetscInt*);
  PetscErrorCode (*VecZeroEntries)(Vec);
  PetscErrorCode (*VecScale)(Vec,PetscScalar);
  PetscErrorCode (*VecSet)(Vec,PetscScalar);
  PetscErrorCode (*VecSetValues)(Vec,PetscInt,const PetscInt[],const PetscScalar[],InsertMode);
  PetscErrorCode (*VecGetValues)(Vec,PetscInt,const PetscInt[],PetscScalar[]);
  PetscErrorCode (*VecAssemblyBegin)(Vec);
  PetscErrorCode (*VecAssemblyEnd)(Vec);
  PetscErrorCode (*VecNorm)(Vec,NormType,PetscReal *);

} petsc_methods_table;

typedef struct
{
  void* petsc;
  petsc_methods_table methods;
  bool finalize_petsc;
} petsc_factory_t;

typedef struct
{
  petsc_factory_t* factory;
  KSP ksp;
} petsc_solver_t;

typedef struct
{
  petsc_factory_t* factory;
  Mat A;
} petsc_matrix_t;

typedef struct
{
  petsc_factory_t* factory;
  Vec v;
} petsc_vector_t;

static void petsc_solver_set_tolerances(void* context,
                                        real_t rel_tol,
                                        real_t abs_tol,
                                        real_t div_tol)

{
  petsc_solver_t* solver = context;
  PetscReal r, a, d;
  PetscInt iters;
  solver->factory->methods.KSPGetTolerances(solver->ksp, &r, &a, &d, &iters);
  solver->factory->methods.KSPSetTolerances(solver->ksp, rel_tol, abs_tol, 
                                            div_tol, iters);
}

static void petsc_solver_set_max_iterations(void* context,
                                            int max_iters)
{
  petsc_solver_t* solver = context;
  PetscReal r, a, d;
  PetscInt iters;
  solver->factory->methods.KSPGetTolerances(solver->ksp, &r, &a, &d, &iters);
  solver->factory->methods.KSPSetTolerances(solver->ksp, r, a, d, max_iters);
}

static void petsc_solver_set_operator(void* context,
                                      void* op)
{
  petsc_solver_t* solver = context;
  Mat A = op;
  PetscErrorCode err = solver->factory->methods.KSPSetOperators(solver->ksp, A, A);
  if (err != 0)
    polymec_error("petsc_solver_set_operator failed!");
}

static bool petsc_solver_solve(void* context,
                               void* x,
                               void* b,
                               real_t* res_norm,
                               int* num_iters)
{
  petsc_solver_t* solver = context;
  Vec X = x;
  Vec B = b;
  PetscErrorCode err = solver->factory->methods.KSPSolve(solver->ksp, X, B);
  if (err != 0)
    polymec_error("petsc_solver_solve: Call to linear solve failed!");
  int code;
  PetscInt niters;
  solver->factory->methods.KSPGetConvergedReason(solver->ksp, &code);
  solver->factory->methods.KSPGetIterationNumber(solver->ksp, &niters);
  *num_iters = (int)niters;
  solver->factory->methods.KSPGetResidualNorm(solver->ksp, res_norm);
  return (code > 0);
}

static void petsc_solver_dtor(void* context)
{
  petsc_solver_t* solver = context;
  solver->factory->methods.KSPDestroy(&solver->ksp);
  solver->factory = NULL;
  polymec_free(solver);
}

static krylov_solver_t* petsc_factory_gmres_solver(void* context,
                                                   MPI_Comm comm,
                                                   int krylov_dimension)
{
  petsc_solver_t* solver = polymec_malloc(sizeof(petsc_solver_t));
  solver->factory = context;
  solver->factory->methods.KSPCreate(comm, &solver->ksp);
  solver->factory->methods.KSPSetType(solver->ksp, "gmres");
  // FIXME: We ignore the Krylov subspace dimension for now.
  // FIXME: Consider altering orthogonalization scheme?

  // Handle the preconditioner's options.
  PC pc;
  solver->factory->methods.KSPGetPC(solver->ksp, &pc);
  solver->factory->methods.PCSetFromOptions(pc);

  // Set the thing up.
  solver->factory->methods.KSPSetFromOptions(solver->ksp);
  solver->factory->methods.KSPSetUp(solver->ksp);

  // Set up the virtual table.
  krylov_solver_vtable vtable = {.set_tolerances = petsc_solver_set_tolerances,
                                 .set_max_iterations = petsc_solver_set_max_iterations,
                                 .set_operator = petsc_solver_set_operator,
                                 .solve = petsc_solver_solve,
                                 .dtor = petsc_solver_dtor};
  return krylov_solver_new("PETSc GMRES", solver, vtable);
}

static krylov_solver_t* petsc_factory_bicgstab_solver(void* context,
                                                      MPI_Comm comm)
{
  petsc_solver_t* solver = polymec_malloc(sizeof(petsc_solver_t));
  solver->factory = context;
  solver->factory->methods.KSPCreate(comm, &solver->ksp);
  solver->factory->methods.KSPSetType(solver->ksp, "bicgstab");

  // Handle the preconditioner's options.
  PC pc;
  solver->factory->methods.KSPGetPC(solver->ksp, &pc);
  solver->factory->methods.PCSetFromOptions(pc);

  // Set the thing up.
  solver->factory->methods.KSPSetFromOptions(solver->ksp);
  solver->factory->methods.KSPSetUp(solver->ksp);

  // Set up the virtual table.
  krylov_solver_vtable vtable = {.set_tolerances = petsc_solver_set_tolerances,
                                 .set_max_iterations = petsc_solver_set_max_iterations,
                                 .set_operator = petsc_solver_set_operator,
                                 .solve = petsc_solver_solve,
                                 .dtor = petsc_solver_dtor};
  return krylov_solver_new("PETSc Bi-CGSTAB", solver, vtable);
}

static krylov_solver_t* petsc_factory_special_solver(void* context,
                                                     MPI_Comm comm,
                                                     const char* solver_name,
                                                     string_string_unordered_map_t* options)
{
  petsc_solver_t* solver = polymec_malloc(sizeof(petsc_solver_t));
  solver->factory = context;
  solver->factory->methods.KSPCreate(comm, &solver->ksp);
  solver->factory->methods.KSPSetType(solver->ksp, solver_name);

  // Handle the preconditioner's options.
  PC pc;
  solver->factory->methods.KSPGetPC(solver->ksp, &pc);
  solver->factory->methods.PCSetFromOptions(pc);

  // Set the thing up.
  solver->factory->methods.KSPSetFromOptions(solver->ksp);
  solver->factory->methods.KSPSetUp(solver->ksp);

  // Set up the virtual table.
  krylov_solver_vtable vtable = {.set_tolerances = petsc_solver_set_tolerances,
                                 .set_max_iterations = petsc_solver_set_max_iterations,
                                 .set_operator = petsc_solver_set_operator,
                                 .solve = petsc_solver_solve,
                                 .dtor = petsc_solver_dtor};
  return krylov_solver_new(solver_name, solver, vtable);
}

static krylov_pc_t* petsc_factory_pc(void* context,
                                     MPI_Comm comm,
                                     const char* pc_name,
                                     string_string_unordered_map_t* options)
{
  petsc_factory_t* factory = context;

  PC pc;
  factory->methods.PCCreate(comm, &pc);
  factory->methods.PCSetType(pc, pc_name);

  // FIXME: Options go here.

  // Set up the virtual table.
  krylov_pc_vtable vtable = {.dtor = NULL}; // FIXME: leak?
  return krylov_pc_new(pc_name, pc, vtable);
}

static void* petsc_matrix_clone(void* context)
{
  petsc_matrix_t* A = context;
  petsc_matrix_t* clone = polymec_malloc(sizeof(petsc_matrix_t));
  clone->factory = A->factory;
  PetscErrorCode err = A->factory->methods.MatConvert(A->A, MATSAME, MAT_INITIAL_MATRIX, &(clone->A));
  if (err != 0)
    polymec_error("petsc_matrix_clone: Error!");
  return clone;
}

static void petsc_matrix_zero(void* context)
{
  petsc_matrix_t* A = context;
  A->factory->methods.MatZeroEntries(A->A);
}

static void petsc_matrix_scale(void* context, real_t scale_factor)
{
  petsc_matrix_t* A = context;
  A->factory->methods.MatScale(A->A, scale_factor);
}

static void petsc_matrix_add_identity(void* context, real_t scale_factor)
{
  petsc_matrix_t* A = context;
  A->factory->methods.MatShift(A->A, scale_factor);
}

static void petsc_matrix_set_diagonal(void* context, void* D)
{
  petsc_matrix_t* A = context;
  Vec diag = D;
  A->factory->methods.MatDiagonalSet(A->A, diag, INSERT_VALUES);
}

static void petsc_matrix_add_diagonal(void* context, void* D)
{
  petsc_matrix_t* A = context;
  Vec diag = D;
  A->factory->methods.MatDiagonalSet(A->A, diag, ADD_VALUES);
}

static void petsc_matrix_enter_values(void* context, index_t num_rows,
                                      index_t* num_columns, index_t* rows, index_t* columns,
                                      real_t* values, InsertMode insert_mode)
{
  petsc_matrix_t* A = context;
  int col_offset = 0;
  for (int r = 0; r < num_rows; ++r)
  {
    int num_cols = num_columns[r];
    PetscInt row = rows[r];
    PetscInt cols[num_cols];
    real_t vals[num_cols];
    for (int c = 0; c < num_cols; ++c, ++col_offset)
    {
      cols[c] = columns[col_offset];
      vals[c] = values[col_offset];
    }
    A->factory->methods.MatSetValues(A->A, 1, &row, num_cols, cols, vals, insert_mode);
  }
}

static void petsc_matrix_set_values(void* context, index_t num_rows,
                                    index_t* num_columns, index_t* rows, index_t* columns,
                                    real_t* values)
{
  petsc_matrix_enter_values(context, num_rows, num_columns, rows, columns, values, INSERT_VALUES);
}

static void petsc_matrix_add_values(void* context, index_t num_rows,
                                    index_t* num_columns, index_t* rows, index_t* columns,
                                    real_t* values)
{
  petsc_matrix_enter_values(context, num_rows, num_columns, rows, columns, values, ADD_VALUES);
}

static void petsc_matrix_start_assembly(void* context)
{
  petsc_matrix_t* A = context;
  A->factory->methods.MatAssemblyBegin(A->A, MAT_FINAL_ASSEMBLY);
}

static void petsc_matrix_finish_assembly(void* context)
{
  petsc_matrix_t* A = context;
  A->factory->methods.MatAssemblyEnd(A->A, MAT_FINAL_ASSEMBLY);
}

static void petsc_matrix_get_values(void* context, index_t num_rows,
                                    index_t* num_columns, index_t* rows, index_t* columns,
                                    real_t* values)
{
  petsc_matrix_t* A = context;
  int col_offset = 0;
  for (int r = 0; r < num_rows; ++r)
  {
    int num_cols = num_columns[r];
    PetscInt row = rows[r];
    PetscInt cols[num_cols];
    for (int c = 0; c < num_cols; ++c, ++col_offset)
      cols[c] = columns[col_offset];
    A->factory->methods.MatGetValues(A->A, 1, &row, num_cols, cols, &values[col_offset-num_cols]);
  }
}

static void petsc_matrix_dtor(void* context)
{
  petsc_matrix_t* A = context;
  A->factory->methods.MatDestroy(&(A->A));
  A->factory = NULL;
  polymec_free(A);
}

static krylov_matrix_t* petsc_factory_matrix(void* context,
                                             adj_graph_t* sparsity)
{
  petsc_matrix_t* A = polymec_malloc(sizeof(petsc_matrix_t));
  A->factory = context;
  int N_local = adj_graph_num_vertices(sparsity);
  index_t* vtx_dist = adj_graph_vertex_dist(sparsity);
  MPI_Comm comm = adj_graph_comm(sparsity);
  int nprocs;
  MPI_Comm_size(comm, &nprocs);
  int N_global = (int)vtx_dist[nprocs];
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

  // Set up the virtual table.
  krylov_matrix_vtable vtable = {.clone = petsc_matrix_clone,
                                 .zero = petsc_matrix_zero,
                                 .scale = petsc_matrix_scale,
                                 .add_identity = petsc_matrix_add_identity,
                                 .add_diagonal = petsc_matrix_add_diagonal,
                                 .set_diagonal = petsc_matrix_set_diagonal,
                                 .set_values = petsc_matrix_set_values,
                                 .add_values = petsc_matrix_add_values,
                                 .start_assembly = petsc_matrix_start_assembly,
                                 .finish_assembly = petsc_matrix_finish_assembly,
                                 .get_values = petsc_matrix_get_values,
                                 .dtor = petsc_matrix_dtor};
  return krylov_matrix_new(A, vtable, N_local, N_global);
}

static krylov_matrix_t* petsc_factory_block_matrix(void* context,
                                                   adj_graph_t* sparsity,
                                                   int block_size)
{
  petsc_matrix_t* A = polymec_malloc(sizeof(petsc_matrix_t));
  A->factory = context;
  int N_local = adj_graph_num_vertices(sparsity);
  index_t* vtx_dist = adj_graph_vertex_dist(sparsity);
  MPI_Comm comm = adj_graph_comm(sparsity);
  int nprocs;
  MPI_Comm_size(comm, &nprocs);
  int N_global = (int)vtx_dist[nprocs];
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

  // Set up the virtual table.
  krylov_matrix_vtable vtable = {.clone = petsc_matrix_clone,
                                 .zero = petsc_matrix_zero,
                                 .scale = petsc_matrix_scale,
                                 .add_identity = petsc_matrix_add_identity,
                                 .add_diagonal = petsc_matrix_add_diagonal,
                                 .set_diagonal = petsc_matrix_set_diagonal,
                                 .set_values = petsc_matrix_set_values,
                                 .add_values = petsc_matrix_add_values,
                                 .start_assembly = petsc_matrix_start_assembly,
                                 .finish_assembly = petsc_matrix_finish_assembly,
                                 .get_values = petsc_matrix_get_values,
                                 .dtor = petsc_matrix_dtor};
  return krylov_matrix_new(A, vtable, N_local, N_global);
}

static void* petsc_vector_clone(void* context)
{
  petsc_vector_t* v = context;
  petsc_vector_t* clone = polymec_malloc(sizeof(petsc_vector_t));
  clone->factory = v->factory;
  PetscErrorCode err = v->factory->methods.VecDuplicate(v->v, &(clone->v));
  if (err != 0)
    polymec_error("petsc_vector_clone: Error in VecDuplicate!");
  err = v->factory->methods.VecCopy(v->v, clone->v);
  if (err != 0)
    polymec_error("petsc_vector_clone: Error in VecCopy!");
  return clone;
}

static void petsc_vector_zero(void* context)
{
  petsc_vector_t* v = context;
  v->factory->methods.VecZeroEntries(v->v);
}

static void petsc_vector_set_value(void* context, real_t value)
{
  petsc_vector_t* v = context;
  v->factory->methods.VecSet(v->v, value);
}

static void petsc_vector_scale(void* context, real_t scale_factor)
{
  petsc_vector_t* v = context;
  v->factory->methods.VecScale(v->v, scale_factor);
}

static void petsc_vector_set_values(void* context, index_t num_values,
                                    index_t* indices, real_t* values)
{
  petsc_vector_t* v = context;
  v->factory->methods.VecSetValues(v->v, num_values, indices, values, INSERT_VALUES);
}

static void petsc_vector_add_values(void* context, index_t num_values,
                                    index_t* indices, real_t* values)
{
  petsc_vector_t* v = context;
  v->factory->methods.VecSetValues(v->v, num_values, indices, values, ADD_VALUES);
}

static void petsc_vector_start_assembly(void* context)
{
  petsc_vector_t* v = context;
  v->factory->methods.VecAssemblyBegin(v->v);
}

static void petsc_vector_finish_assembly(void* context)
{
  petsc_vector_t* v = context;
  v->factory->methods.VecAssemblyEnd(v->v);
}

static void petsc_vector_get_values(void* context, index_t num_values,
                                    index_t* indices, real_t* values)
{
  petsc_vector_t* v = context;
  v->factory->methods.VecGetValues(v->v, num_values, indices, values);
}

static real_t petsc_vector_norm(void* context, int p)
{
  petsc_vector_t* v = context;
  real_t norm;
  NormType norm_type;
  if (p == 0)
    norm_type = NORM_INFINITY;
  else if (p == 1)
    norm_type = NORM_1;
  else
    norm_type = NORM_2;
  v->factory->methods.VecNorm(v->v, norm_type, &norm);
  return norm;
}

static void petsc_vector_dtor(void* context)
{
  petsc_vector_t* v = context;
  v->factory->methods.VecDestroy(&(v->v));
  v->factory = NULL;
  polymec_free(v);
}

static krylov_vector_t* petsc_factory_vector(void* context,
                                             adj_graph_t* dist_graph)
{
  petsc_vector_t* v = polymec_malloc(sizeof(petsc_vector_t));
  v->factory = context;
  MPI_Comm comm = adj_graph_comm(dist_graph);
  int rank, nprocs;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &nprocs);
  index_t* vtx_dist = adj_graph_vertex_dist(dist_graph);
  if (nprocs == 1)
    v->factory->methods.VecCreateSeq(comm, vtx_dist[1], &v->v);
  else
    v->factory->methods.VecCreateMPI(comm, vtx_dist[nprocs], vtx_dist[rank+1]-vtx_dist[rank], &v->v);

  // Set up the virtual table.
  krylov_vector_vtable vtable = {.clone = petsc_vector_clone,
                                 .zero = petsc_vector_zero,
                                 .set_value = petsc_vector_set_value,
                                 .scale = petsc_vector_scale,
                                 .set_values = petsc_vector_set_values,
                                 .add_values = petsc_vector_add_values,
                                 .start_assembly = petsc_vector_start_assembly,
                                 .finish_assembly = petsc_vector_finish_assembly,
                                 .get_values = petsc_vector_get_values,
                                 .norm = petsc_vector_norm,
                                 .dtor = petsc_vector_dtor};
  return krylov_vector_new(v, vtable, vtx_dist[rank+1]-vtx_dist[rank], vtx_dist[nprocs]);
}

static void petsc_factory_dtor(void* context)
{
  petsc_factory_t* factory = context;
  if (factory->finalize_petsc)
  {
    log_debug("petsc_krylov_factory: Finalizing PETSc.");
    factory->methods.PetscFinalize();
  }
  log_debug("petsc_krylov_factory: Closing PETSc library.");
  dlclose(factory->petsc);
  polymec_free(factory);
}

// Use this to retrieve symbols from dynamically loaded libraries.
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

krylov_factory_t* petsc_krylov_factory(const char* petsc_dir,
                                       const char* petsc_arch)
{
  ASSERT(((petsc_dir == NULL) && (petsc_arch == NULL)) ||
         ((petsc_dir != NULL) && (petsc_arch != NULL)));

  petsc_factory_t* factory = polymec_malloc(sizeof(petsc_factory_t));

  // Try to find PETSc.
  char petsc_path[FILENAME_MAX+1];
  char *my_petsc_dir, *my_petsc_arch;
  if (petsc_dir == NULL)
  {
    my_petsc_dir = getenv("PETSC_DIR");
    my_petsc_arch = getenv("PETSC_ARCH");
  }
  if ((petsc_dir == NULL) && (my_petsc_dir == NULL))
    polymec_error("PETSC directory was not given. Please set PETSC_DIR.");
  if ((petsc_arch == NULL) && (my_petsc_arch == NULL))
    polymec_error("PETSC architecture was not given. Please set PETSC_ARCH.");
  if (petsc_dir != NULL)
    snprintf(petsc_path, FILENAME_MAX, "%s/%s/lib/libpetsc%s", petsc_dir, petsc_arch, SHARED_LIBRARY_SUFFIX);
  else
    snprintf(petsc_path, FILENAME_MAX, "%s/%s/lib/libpetsc%s", my_petsc_dir, my_petsc_arch, SHARED_LIBRARY_SUFFIX);

  // Try to open libPETSc and mine it for symbols.
  log_debug("petsc_krylov_factory: Opening PETSc library at %s.", petsc_path);
  void* petsc = dlopen(petsc_path, RTLD_NOW);
  if (petsc == NULL)
  {
    char* msg = dlerror();
    polymec_error("petsc_krylov_factory: %s.", msg);
  }

  // Fetch the version.
  char version[128];
  PetscErrorCode (*PetscGetVersion)(char* version, size_t len);
  FETCH_SYMBOL(petsc, "PetscGetVersion", PetscGetVersion, failure);
  PetscErrorCode err = PetscGetVersion(version, 128);
  if (err != 0)
  {
    log_urgent("petsc_krylov_factory: Error getting PETSc version string.");
    goto failure;
  }

  log_debug("petsc_krylov_factory: Got PETSc version: %s", version);

  // In the absence of another mechanism to determine whether PETSc was 
  // configured to use 64-bit indexing, we read through petscconf.h
  // to check for the presence of PETSC_USE_64BIT_INDICES.
  bool petsc_uses_64bit_indices = false;
  {
    char petsc_arch_incdir[FILENAME_MAX];
    snprintf(petsc_arch_incdir, FILENAME_MAX, "%s/%s/include", petsc_dir, petsc_arch);
    if (!directory_exists(petsc_arch_incdir))
    {
      log_urgent("petsc_krylov_factory: Could not find directory %s.", petsc_arch_incdir);
      goto failure;
    }
    char petscconf_h[FILENAME_MAX];
    snprintf(petscconf_h, FILENAME_MAX, "%s/petscconf.h", petsc_arch_incdir);
    text_buffer_t* buffer = text_buffer_from_file(petscconf_h);
    if (buffer == NULL)
    {
      log_urgent("petsc_krylov_factory: Could not read petscconf.h.");
      goto failure;
    }

    // Look for PETSC_USE_64BIT_INDICES.
    int pos = 0, length;
    char* line;
    while (text_buffer_next_nonempty(buffer, &pos, &line, &length))
    {
      char text[length+1];
      string_copy_from_raw(line, length, text);
      if (string_contains(text, "PETSC_USE_64BIT_INDICES"))
      {
        petsc_uses_64bit_indices = true;
        break;
      }
    }
    text_buffer_free(buffer);
  }
#if 0
  if (sizeof(index_t) == sizeof(int64_t))
  {
    if (!petsc_uses_64bit_indices)
    {
      log_urgent("petsc_krylov_factory: Since polymec is configured for 64-bit indices,\n"
                 "  PETSc must be built using --with-64-bit-indices.");
      goto failure;
    }
  }
  else
  {
    if (petsc_uses_64bit_indices)
    {
      log_urgent("petsc_krylov_factory: Since polymec is configured for 32-bit indices,\n"
                 "  PETSc must be built with 32-bit indices (not using --with-64-bit-indices).");
      goto failure;
    }
  }
#endif
  if (petsc_uses_64bit_indices)
  {
    log_urgent("petsc_krylov_factory: Currently, PETSc must be configured to use 32-bit indices.");
    goto failure;
  }

  // Get the symbols.
#define FETCH_PETSC_SYMBOL(symbol_name) \
  FETCH_SYMBOL(petsc, #symbol_name, factory->methods.symbol_name, failure);

  FETCH_PETSC_SYMBOL(PetscInitialize);
  FETCH_PETSC_SYMBOL(PetscInitialized);
  FETCH_PETSC_SYMBOL(PetscFinalize);

  FETCH_PETSC_SYMBOL(KSPCreate);
  FETCH_PETSC_SYMBOL(KSPSetType);
  FETCH_PETSC_SYMBOL(KSPSetFromOptions);
  FETCH_PETSC_SYMBOL(KSPSetUp);
  FETCH_PETSC_SYMBOL(KSPGetPC);
  FETCH_PETSC_SYMBOL(KSPSetTolerances);
  FETCH_PETSC_SYMBOL(KSPGetTolerances);
  FETCH_PETSC_SYMBOL(KSPSetOperators);
  FETCH_PETSC_SYMBOL(KSPSolve);
  FETCH_PETSC_SYMBOL(KSPGetConvergedReason);
  FETCH_PETSC_SYMBOL(KSPGetIterationNumber);
  FETCH_PETSC_SYMBOL(KSPGetResidualNorm);
  FETCH_PETSC_SYMBOL(KSPDestroy);

  FETCH_PETSC_SYMBOL(PCCreate);
  FETCH_PETSC_SYMBOL(PCSetType);
  FETCH_PETSC_SYMBOL(PCSetFromOptions);

  FETCH_PETSC_SYMBOL(MatCreate);
  FETCH_PETSC_SYMBOL(MatConvert);
  FETCH_PETSC_SYMBOL(MatSetType);
  FETCH_PETSC_SYMBOL(MatSetSizes);
  FETCH_PETSC_SYMBOL(MatSeqAIJSetPreallocation);
  FETCH_PETSC_SYMBOL(MatMPIAIJSetPreallocation);
  FETCH_PETSC_SYMBOL(MatSeqBAIJSetPreallocation);
  FETCH_PETSC_SYMBOL(MatMPIBAIJSetPreallocation);
  FETCH_PETSC_SYMBOL(MatSetBlockSize);
  FETCH_PETSC_SYMBOL(MatSetUp);
  FETCH_PETSC_SYMBOL(MatDestroy);
  FETCH_PETSC_SYMBOL(MatScale);
  FETCH_PETSC_SYMBOL(MatShift);
  FETCH_PETSC_SYMBOL(MatDiagonalSet);
  FETCH_PETSC_SYMBOL(MatZeroEntries);
  FETCH_PETSC_SYMBOL(MatGetSize);
  FETCH_PETSC_SYMBOL(MatGetLocalSize);
  FETCH_PETSC_SYMBOL(MatSetValues);
  FETCH_PETSC_SYMBOL(MatGetValues);
  FETCH_PETSC_SYMBOL(MatAssemblyBegin);
  FETCH_PETSC_SYMBOL(MatAssemblyEnd);
  
  FETCH_PETSC_SYMBOL(VecCreateSeq);
  FETCH_PETSC_SYMBOL(VecCreateMPI);
  FETCH_PETSC_SYMBOL(VecDuplicate);
  FETCH_PETSC_SYMBOL(VecCopy);
  FETCH_PETSC_SYMBOL(VecSetUp);
  FETCH_PETSC_SYMBOL(VecDestroy);
  FETCH_PETSC_SYMBOL(VecGetSize);
  FETCH_PETSC_SYMBOL(VecZeroEntries);
  FETCH_PETSC_SYMBOL(VecScale);
  FETCH_PETSC_SYMBOL(VecSet);
  FETCH_PETSC_SYMBOL(VecSetValues);
  FETCH_PETSC_SYMBOL(VecGetValues);
  FETCH_PETSC_SYMBOL(VecAssemblyBegin);
  FETCH_PETSC_SYMBOL(VecAssemblyEnd);
  FETCH_PETSC_SYMBOL(VecNorm);
#undef FETCH_PETSC_SYMBOL

  log_debug("petsc_krylov_factory: Got PETSc symbols.");

  // Initialize PETSc if needed.
  PetscBool initialized;
  factory->methods.PetscInitialized(&initialized);
  if (initialized == PETSC_FALSE)
  {
    log_debug("petsc_krylov_factory: Initializing PETSc.");
    // Build a list of command line arguments.
    options_t* opts = options_argv();
    int argc = options_num_arguments(opts);
    char** argv = polymec_malloc(sizeof(char*));
    for (int i = 0; i < argc; ++i)
      argv[i] = options_argument(opts, i);
    factory->methods.PetscInitialize(&argc, &argv, NULL, NULL);
    polymec_free(argv);
    factory->finalize_petsc = true;
  }

  // Stash the library.
  factory->petsc = petsc; 

  // Construct the factory.
  krylov_factory_vtable vtable = {.gmres_solver = petsc_factory_gmres_solver,
                                  .bicgstab_solver = petsc_factory_bicgstab_solver,
                                  .special_solver = petsc_factory_special_solver,
                                  .preconditioner = petsc_factory_pc,
                                  .matrix = petsc_factory_matrix,
                                  .block_matrix = petsc_factory_block_matrix,
                                  .vector = petsc_factory_vector,
                                  .dtor = petsc_factory_dtor};
  return krylov_factory_new(version, factory, vtable);

failure:
  dlclose(petsc);
  polymec_free(factory);
  return NULL;
}

