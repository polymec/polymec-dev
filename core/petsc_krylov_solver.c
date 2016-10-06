// Copyright (c) 2012-2016, Jeffrey N. Johnson
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

#if POLYMEC_HAVE_SHARED_LIBS && POLYMEC_HAVE_MPI

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
typedef void* PetscViewer;
typedef enum {NORM_1=0,NORM_2=1,NORM_FROBENIUS=2,NORM_INFINITY=3,NORM_1_AND_2=4} NormType;
typedef enum {KSP_NORM_DEFAULT = -1,KSP_NORM_NONE = 0,KSP_NORM_PRECONDITIONED = 1,KSP_NORM_UNPRECONDITIONED = 2,KSP_NORM_NATURAL = 3} KSPNormType;
typedef enum {MAT_FLUSH_ASSEMBLY=1,MAT_FINAL_ASSEMBLY=0} MatAssemblyType;
typedef enum {INSERT_VALUES=1,ADD_VALUES=0} InsertMode;
typedef enum {MAT_INITIAL_MATRIX,MAT_REUSE_MATRIX,MAT_IGNORE_MATRIX} MatReuse;
typedef enum {DIFFERENT_NONZERO_PATTERN,SUBSET_NONZERO_PATTERN,SAME_NONZERO_PATTERN} MatStructure;
static const char* MATSEQAIJ = "seqaij";
static const char* MATMPIAIJ = "mpiaij";
static const char* MATSEQBAIJ = "seqbaij";
static const char* MATMPIBAIJ = "mpibaij";
static const char* MATSAME = "same";
#define PETSC_DEFAULT -2
#define PETSC_DECIDE -1
#define PETSC_DETERMINE PETSC_DECIDE
typedef struct
{
  PetscErrorCode (*PetscInitializeNoArguments)();
  PetscErrorCode (*PetscInitialized)(PetscBool*);
  PetscErrorCode (*PetscFinalize)();
  PetscErrorCode (*PetscPushErrorHandler)(PetscErrorCode (*handler)(MPI_Comm comm, int, const char*, const char*, PetscErrorCode, int, const char*, void*), void*);
  void* (*PETSC_VIEWER_STDOUT_)(MPI_Comm);
  PetscErrorCode (*PetscViewerASCIIOpenWithFILE)(MPI_Comm, FILE*, PetscViewer*);
  PetscErrorCode (*PetscViewerDestroy)(PetscViewer*);

  PetscErrorCode (*KSPCreate)(MPI_Comm,KSP *);
  PetscErrorCode (*KSPSetType)(KSP,KSPType);
  PetscErrorCode (*KSPGetPC)(KSP,PC*);
  PetscErrorCode (*KSPSetPC)(KSP,PC);
  PetscErrorCode (*KSPSetFromOptions)(KSP);
  PetscErrorCode (*KSPSetTolerances)(KSP,PetscReal,PetscReal,PetscReal,PetscInt);
  PetscErrorCode (*KSPGetTolerances)(KSP,PetscReal*,PetscReal*,PetscReal*,PetscInt*);
  PetscErrorCode (*KSPSetDiagonalScale)(KSP,PetscBool);
  PetscErrorCode (*KSPSetDiagonalScaleFix)(KSP,PetscBool);
  PetscErrorCode (*KSPSetNormType)(KSP,KSPNormType);
  PetscErrorCode (*KSPSetOperators)(KSP,Mat,Mat);
  PetscErrorCode (*KSPSetUp)(KSP);
  PetscErrorCode (*KSPSolve)(KSP,Vec,Vec);
  PetscErrorCode (*KSPGetConvergedReason)(KSP,int*);
  PetscErrorCode (*KSPGetIterationNumber)(KSP,PetscInt*);
  PetscErrorCode (*KSPGetResidualNorm)(KSP,PetscReal*);
  PetscErrorCode (*KSPDestroy)(KSP*);
  PetscErrorCode (*KSPGMRESSetRestart)(KSP,PetscInt);

  PetscErrorCode (*PCCreate)(MPI_Comm,PC*);
  PetscErrorCode (*PCSetType)(PC, const char*);
  PetscErrorCode (*PCSetFromOptions)(PC);

  PetscErrorCode (*MatCreate)(MPI_Comm,Mat*);
  PetscErrorCode (*MatConvert)(Mat, MatType, MatReuse, Mat*);
  PetscErrorCode (*MatCopy)(Mat, Mat, MatStructure);
  PetscErrorCode (*MatSetType)(Mat, MatType);
  PetscErrorCode (*MatSetSizes)(Mat, PetscInt m, PetscInt n, PetscInt M, PetscInt N);
  PetscErrorCode (*MatGetOwnershipRange)(Mat, PetscInt* m, PetscInt* n);
  PetscErrorCode (*MatSeqAIJSetPreallocation)(Mat, PetscInt, PetscInt[]);
  PetscErrorCode (*MatMPIAIJSetPreallocation)(Mat, PetscInt, PetscInt[], PetscInt, PetscInt[]);
  PetscErrorCode (*MatSeqBAIJSetPreallocation)(Mat, PetscInt, PetscInt, PetscInt[]);
  PetscErrorCode (*MatMPIBAIJSetPreallocation)(Mat, PetscInt, PetscInt, PetscInt[], PetscInt, PetscInt[]);
  PetscErrorCode (*MatSetBlockSize)(Mat, PetscInt);
  PetscErrorCode (*MatDestroy)(Mat*);
  PetscErrorCode (*MatScale)(Mat,PetscScalar);
  PetscErrorCode (*MatDiagonalScale)(Mat,Vec,Vec);
  PetscErrorCode (*MatShift)(Mat,PetscScalar);
  PetscErrorCode (*MatDiagonalSet)(Mat,Vec,InsertMode);
  PetscErrorCode (*MatZeroEntries)(Mat);
  PetscErrorCode (*MatMult)(Mat,Vec,Vec);
  PetscErrorCode (*MatMultTranspose)(Mat,Vec,Vec);
  PetscErrorCode (*MatGetSize)(Mat,PetscInt*,PetscInt*);
  PetscErrorCode (*MatGetLocalSize)(Mat,PetscInt*,PetscInt*);
  PetscErrorCode (*MatSetValues)(Mat,PetscInt,const PetscInt[],PetscInt,const PetscInt[],const PetscScalar[],InsertMode);
  PetscErrorCode (*MatAssemblyBegin)(Mat,MatAssemblyType);
  PetscErrorCode (*MatAssemblyEnd)(Mat,MatAssemblyType);
  PetscErrorCode (*MatGetValues)(Mat,PetscInt,const PetscInt[],PetscInt,const PetscInt[],PetscScalar[]);
  PetscErrorCode (*MatView)(Mat,PetscViewer);

  PetscErrorCode (*VecCreateSeq)(MPI_Comm,PetscInt,Vec*);
  PetscErrorCode (*VecCreateMPI)(MPI_Comm,PetscInt,PetscInt,Vec*);
  PetscErrorCode (*VecDuplicate)(Vec,Vec*);
  PetscErrorCode (*VecCopy)(Vec,Vec);
  PetscErrorCode (*VecSetUp)(Vec);
  PetscErrorCode (*VecDestroy)(Vec*);
  PetscErrorCode (*VecGetSize)(Vec, PetscInt*);
  PetscErrorCode (*VecGetLocalSize)(Vec, PetscInt*);
  PetscErrorCode (*VecZeroEntries)(Vec);
  PetscErrorCode (*VecScale)(Vec,PetscScalar);
  PetscErrorCode (*VecSet)(Vec,PetscScalar);
  PetscErrorCode (*VecSetValues)(Vec,PetscInt,const PetscInt[],const PetscScalar[],InsertMode);
  PetscErrorCode (*VecGetValues)(Vec,PetscInt,const PetscInt[],PetscScalar[]);
  PetscErrorCode (*VecGetArray)(Vec,PetscScalar**);
  PetscErrorCode (*VecRestoreArray)(Vec,PetscScalar**);
  PetscErrorCode (*VecAssemblyBegin)(Vec);
  PetscErrorCode (*VecAssemblyEnd)(Vec);
  PetscErrorCode (*VecDot)(Vec,Vec,PetscReal *);
  PetscErrorCode (*VecNorm)(Vec,NormType,PetscReal *);
  PetscErrorCode (*VecView)(Vec,PetscViewer);

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
  MPI_Comm comm;
  Mat A;
  int block_size;
  int* block_sizes;
} petsc_matrix_t;

typedef struct
{
  petsc_factory_t* factory;
  MPI_Comm comm;
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

static inline void petsc_matrix_assemble(void* context)
{
  petsc_matrix_t* A = context;
  A->factory->methods.MatAssemblyBegin(A->A, MAT_FINAL_ASSEMBLY);
  A->factory->methods.MatAssemblyEnd(A->A, MAT_FINAL_ASSEMBLY);
}

static inline void petsc_vector_assemble(void* context)
{
  petsc_vector_t* v = context;
  v->factory->methods.VecAssemblyBegin(v->v);
  v->factory->methods.VecAssemblyEnd(v->v);
}

static void petsc_solver_set_operator(void* context,
                                      void* op)
{
  petsc_solver_t* solver = context;
  petsc_matrix_t* A = op;
  PetscErrorCode err = solver->factory->methods.KSPSetOperators(solver->ksp, A->A, A->A);
  if (err != 0)
    polymec_error("petsc_solver_set_operator failed setting operator.");
  err = solver->factory->methods.KSPSetUp(solver->ksp);
  if (err != 0)
    polymec_error("petsc_solver_set_operator failed setting up solver.");
}

static void petsc_solver_set_pc(void* context,
                                void* pc)
{
  petsc_solver_t* solver = context;
  PC p = pc;
  PetscErrorCode err = solver->factory->methods.KSPSetPC(solver->ksp, p);
  if (err != 0)
    polymec_error("petsc_solver_set_pc failed!");
}

static bool petsc_solver_solve(void* context,
                               void* b,
                               void* x,
                               real_t* res_norm,
                               int* num_iters)
{
  petsc_solver_t* solver = context;
  petsc_vector_t* B = b;
  petsc_vector_t* X = x;
  PetscErrorCode err = solver->factory->methods.KSPSolve(solver->ksp, B->v, X->v);
  if (err != 0)
    polymec_error("petsc_solver_solve: Call to linear solve failed!");
  int code;
  PetscInt niters;
  solver->factory->methods.KSPGetConvergedReason(solver->ksp, &code);
  solver->factory->methods.KSPGetIterationNumber(solver->ksp, &niters);
  *num_iters = (int)niters;
  solver->factory->methods.KSPGetResidualNorm(solver->ksp, res_norm);
  if ((code < 0) && (log_level() == LOG_DEBUG))
  {
    char reason[256]; 
    switch(code)
    {
      case -2: sprintf(reason, "KSP_DIVERGED_NULL"); break;
      case -3: sprintf(reason, "max iterations exceeded"); break;
      case -4: sprintf(reason, "residual norm exceeds divergence tolerance"); break;
      case -5: sprintf(reason, "Krylov method breakdown (singular A or P?)"); break;
      case -6: sprintf(reason, "Krylov method breakdown in biconjugate gradient method"); break;
      case -7: sprintf(reason, "Nonsymmetric matrix given for symmetric Krylov method"); break;
      case -8: sprintf(reason, "Indefinite preconditioner for method requiring SPD matrix"); break;
      case -9: sprintf(reason, "NaN or Inf value encountered in linear solve"); break;
      case -10: sprintf(reason, "Indefinite matrix given for method requiring SPD matrix"); break;
      case -11: sprintf(reason, "Setup for preconditioner failed"); break;
      default: sprintf(reason, "reason unknown");
    }
    log_debug("petsc_solver_solve: Linear solve failed: %s", reason);
  }
  return (code > 0);
}

static void petsc_solver_dtor(void* context)
{
  petsc_solver_t* solver = context;
  solver->factory->methods.KSPDestroy(&solver->ksp);
  solver->factory = NULL;
  polymec_free(solver);
}

static krylov_solver_t* petsc_factory_pcg_solver(void* context,
                                                 MPI_Comm comm)
{
  petsc_solver_t* solver = polymec_malloc(sizeof(petsc_solver_t));
  solver->factory = context;
  solver->factory->methods.KSPCreate(comm, &solver->ksp);
  solver->factory->methods.KSPSetType(solver->ksp, "cg");

  // We always use the unpreconditioned residual norm to determine 
  // convergence.
  solver->factory->methods.KSPSetNormType(solver->ksp, KSP_NORM_UNPRECONDITIONED);

  // Handle the preconditioner's options.
  PC pc;
  solver->factory->methods.KSPGetPC(solver->ksp, &pc);
  solver->factory->methods.PCSetFromOptions(pc);

  // Set the thing up.
  solver->factory->methods.KSPSetFromOptions(solver->ksp);

  // Set up the virtual table.
  krylov_solver_vtable vtable = {.set_tolerances = petsc_solver_set_tolerances,
                                 .set_max_iterations = petsc_solver_set_max_iterations,
                                 .set_operator = petsc_solver_set_operator,
                                 .set_preconditioner = petsc_solver_set_pc,
                                 .solve = petsc_solver_solve,
                                 .dtor = petsc_solver_dtor};
  return krylov_solver_new("PETSc PCG", solver, vtable);
}

static krylov_solver_t* petsc_factory_gmres_solver(void* context,
                                                   MPI_Comm comm,
                                                   int krylov_dimension)
{
  petsc_solver_t* solver = polymec_malloc(sizeof(petsc_solver_t));
  solver->factory = context;
  solver->factory->methods.KSPCreate(comm, &solver->ksp);
  solver->factory->methods.KSPSetType(solver->ksp, "gmres");

  // We always use the unpreconditioned residual norm to determine 
  // convergence.
  solver->factory->methods.KSPSetNormType(solver->ksp, KSP_NORM_UNPRECONDITIONED);

  solver->factory->methods.KSPGMRESSetRestart(solver->ksp, (PetscInt)krylov_dimension);
  // FIXME: Consider altering orthogonalization scheme?

  // Handle the preconditioner's options.
  PC pc;
  solver->factory->methods.KSPGetPC(solver->ksp, &pc);
  solver->factory->methods.PCSetFromOptions(pc);

  // Set the thing up.
  solver->factory->methods.KSPSetFromOptions(solver->ksp);

  // Set up the virtual table.
  krylov_solver_vtable vtable = {.set_tolerances = petsc_solver_set_tolerances,
                                 .set_max_iterations = petsc_solver_set_max_iterations,
                                 .set_operator = petsc_solver_set_operator,
                                 .set_preconditioner = petsc_solver_set_pc,
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
  solver->factory->methods.KSPSetType(solver->ksp, "bcgs");

  // We always use the unpreconditioned residual norm to determine 
  // convergence.
  solver->factory->methods.KSPSetNormType(solver->ksp, KSP_NORM_UNPRECONDITIONED);

  // Handle the preconditioner's options.
  PC pc;
  solver->factory->methods.KSPGetPC(solver->ksp, &pc);
  solver->factory->methods.PCSetFromOptions(pc);

  // Set the thing up.
  solver->factory->methods.KSPSetFromOptions(solver->ksp);

  // Set up the virtual table.
  krylov_solver_vtable vtable = {.set_tolerances = petsc_solver_set_tolerances,
                                 .set_max_iterations = petsc_solver_set_max_iterations,
                                 .set_operator = petsc_solver_set_operator,
                                 .set_preconditioner = petsc_solver_set_pc,
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

  // We always use the unpreconditioned residual norm to determine 
  // convergence.
  solver->factory->methods.KSPSetNormType(solver->ksp, KSP_NORM_UNPRECONDITIONED);

  // Handle the preconditioner's options.
  PC pc;
  solver->factory->methods.KSPGetPC(solver->ksp, &pc);
  solver->factory->methods.PCSetFromOptions(pc);

  // Set the thing up.
  solver->factory->methods.KSPSetFromOptions(solver->ksp);

  // Set up the virtual table.
  krylov_solver_vtable vtable = {.set_tolerances = petsc_solver_set_tolerances,
                                 .set_max_iterations = petsc_solver_set_max_iterations,
                                 .set_operator = petsc_solver_set_operator,
                                 .set_preconditioner = petsc_solver_set_pc,
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
  clone->comm = A->comm;
  clone->block_size = A->block_size;
  if (A->block_sizes != NULL)
  {
    PetscInt N_local, N_local_col;
    A->factory->methods.MatGetLocalSize(A->A, &N_local, &N_local_col);
    clone->block_sizes = polymec_malloc(sizeof(int) * N_local);
    memcpy(clone->block_sizes, A->block_sizes, sizeof(int) * N_local);
  }
  else
    clone->block_sizes = NULL;
  PetscErrorCode err = A->factory->methods.MatConvert(A->A, MATSAME, MAT_INITIAL_MATRIX, &(clone->A));
  if (err != 0)
    polymec_error("petsc_matrix_clone: Error!");
  return clone;
}

static void petsc_matrix_copy(void* context, void* copy)
{
  petsc_matrix_t* A = context;
  petsc_matrix_t* B = copy;
  PetscErrorCode err = A->factory->methods.MatCopy(A->A, B->A, SAME_NONZERO_PATTERN);
  if (err != 0)
    polymec_error("petsc_matrix_copy: Error!");
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

static void petsc_matrix_diag_scale(void* context, 
                                    void* L,
                                    void* R)
{
  petsc_matrix_t* A = context;
  petsc_vector_t* LL = L;
  petsc_vector_t* RR = R;
  void* l = (LL == NULL) ? NULL : LL->v;
  void* r = (RR == NULL) ? NULL : RR->v;
  A->factory->methods.MatDiagonalScale(A->A, l, r);
}

static void petsc_matrix_add_identity(void* context, real_t scale_factor)
{
  petsc_matrix_t* A = context;
  A->factory->methods.MatShift(A->A, scale_factor);
}

static void petsc_matrix_set_diagonal(void* context, void* D)
{
  petsc_matrix_t* A = context;
  petsc_vector_t* diag = D;
  A->factory->methods.MatDiagonalSet(A->A, diag->v, INSERT_VALUES);
}

static void petsc_matrix_add_diagonal(void* context, void* D)
{
  petsc_matrix_t* A = context;
  petsc_vector_t* diag = D;
  A->factory->methods.MatDiagonalSet(A->A, diag->v, ADD_VALUES);
}

static void petsc_matrix_matvec(void* context, void* X, bool transpose, void* Y)
{
  petsc_matrix_t* A = context;
  petsc_vector_t* x = X;
  petsc_vector_t* y = Y;
  if (transpose)
    A->factory->methods.MatMultTranspose(A->A, x->v, y->v);
  else
    A->factory->methods.MatMult(A->A, x->v, y->v);
}

static void petsc_matrix_insert_values(void* context, index_t num_rows,
                                       index_t* num_columns, index_t* rows, index_t* columns,
                                       real_t* values, InsertMode insert_mode)
{
  petsc_matrix_t* A = context;
  int col_offset = 0;
  for (int r = 0; r < num_rows; ++r)
  {
    PetscInt num_cols = num_columns[r];
    PetscInt row = rows[r];
    PetscInt cols[num_cols];
    PetscReal vals[num_cols];
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
  petsc_matrix_insert_values(context, num_rows, num_columns, rows, columns, values, INSERT_VALUES);
}

static void petsc_matrix_add_values(void* context, index_t num_rows,
                                    index_t* num_columns, index_t* rows, index_t* columns,
                                    real_t* values)
{
  petsc_matrix_insert_values(context, num_rows, num_columns, rows, columns, values, ADD_VALUES);
}

static void petsc_matrix_get_values(void* context, index_t num_rows,
                                    index_t* num_columns, index_t* rows, index_t* columns,
                                    real_t* values)
{
  petsc_matrix_t* A = context;
  int col_offset = 0;
  for (int r = 0; r < num_rows; ++r)
  {
    PetscInt num_cols = num_columns[r];
    PetscInt row = rows[r];
    PetscInt cols[num_cols];
    for (int c = 0; c < num_cols; ++c, ++col_offset)
      cols[c] = columns[col_offset];
    PetscReal vals[num_cols];
    A->factory->methods.MatGetValues(A->A, 1, &row, num_cols, cols, vals);
    for (int c = 0; c < num_cols; ++c)
      values[col_offset-num_cols+c] = vals[c];
  }
}

static void petsc_matrix_fprintf(void* context, FILE* stream)
{
  petsc_matrix_t* A = context;
  PetscViewer viewer;
  A->factory->methods.PetscViewerASCIIOpenWithFILE(A->comm, stream, &viewer);
  A->factory->methods.MatView(A->A, viewer);
  A->factory->methods.PetscViewerDestroy(&viewer);
}

static void petsc_matrix_dtor(void* context)
{
  petsc_matrix_t* A = context;
  if (A->block_sizes != NULL)
    polymec_free(A->block_sizes);
  A->factory->methods.MatDestroy(&(A->A));
  A->factory = NULL;
  polymec_free(A);
}

static krylov_matrix_t* petsc_factory_matrix(void* context,
                                             matrix_sparsity_t* sparsity)
{
  petsc_matrix_t* A = polymec_malloc(sizeof(petsc_matrix_t));
  A->factory = context;
  A->comm = matrix_sparsity_comm(sparsity);
  A->block_size = 1;
  A->block_sizes = NULL;

  index_t N_local = matrix_sparsity_num_local_rows(sparsity);
  index_t N_global = matrix_sparsity_num_global_rows(sparsity);

  int nprocs, rank;
  MPI_Comm_size(A->comm, &nprocs);
  MPI_Comm_rank(A->comm, &rank);
  A->factory->methods.MatCreate(A->comm, &A->A);
  if (nprocs == 1)
    A->factory->methods.MatSetType(A->A, MATSEQAIJ);
  else
    A->factory->methods.MatSetType(A->A, MATMPIAIJ);
  A->factory->methods.MatSetSizes(A->A, N_local, N_local, N_global, N_global);

  index_t* row_dist = matrix_sparsity_row_distribution(sparsity);
  PetscInt start = row_dist[rank], end = row_dist[rank+1];

  // Perform preallocation.
  if (nprocs == 1)
  {
    PetscInt nnz[N_local];
    int pos = 0;
    index_t row, r = 0;
    while (matrix_sparsity_next_row(sparsity, &pos, &row))
    {
      nnz[r] = (PetscInt)(matrix_sparsity_num_columns(sparsity, row));
      ++r;
    }
    A->factory->methods.MatSeqAIJSetPreallocation(A->A, PETSC_DEFAULT, nnz);
  }
  else
  {
    // d_nnz[r] is the number of nonzeros for the "diagonal" portion of row r, 
    // and o_nnz[r] is the number of nonzeros for the "off-diagonal" portion.
    PetscInt d_nnz[N_local], o_nnz[N_local];
    memset(d_nnz, 0, sizeof(PetscInt) * N_local);
    memset(o_nnz, 0, sizeof(PetscInt) * N_local);
    index_t row, r = 0;
    int rpos = 0;
    while (matrix_sparsity_next_row(sparsity, &rpos, &row))
    {
      int cpos = 0;
      index_t column;
      while (matrix_sparsity_next_column(sparsity, row, &cpos, &column))
      {
        if ((column < start) || (column >= end))
          o_nnz[r] += 1;
        else
          d_nnz[r] += 1;
      }
      ++r;
    }
    A->factory->methods.MatMPIAIJSetPreallocation(A->A, PETSC_DEFAULT, d_nnz, PETSC_DEFAULT, o_nnz);
  }

  // Set the "non-zero" values to zero initially. This constructs the specific non-zero structure.
  index_t num_columns[N_local], rows[N_local];
  index_t row, r = 0;
  int rpos = 0;
  while (matrix_sparsity_next_row(sparsity, &rpos, &row))
  {
    rows[r] = row;
    num_columns[r] = matrix_sparsity_num_columns(sparsity, row);
    ++r;
  }

  index_t nnz = matrix_sparsity_num_nonzeros(sparsity);
  index_t col_offset = 0;
  index_t columns[nnz];
  rpos = 0;
  while (matrix_sparsity_next_row(sparsity, &rpos, &row))
  {
    int cpos = 0;
    index_t column;
    while (matrix_sparsity_next_column(sparsity, row, &cpos, &column))
    {
      columns[col_offset] = column;
      ++col_offset;
    }
  }
  ASSERT(col_offset == nnz);

  real_t zeros[nnz];
  memset(zeros, 0, sizeof(real_t) * nnz);
  petsc_matrix_set_values(A, N_local, num_columns, rows, columns, zeros);

  // Assemble the matrix.
  petsc_matrix_assemble(A);

  // Set up the virtual table.
  krylov_matrix_vtable vtable = {.clone = petsc_matrix_clone,
                                 .copy = petsc_matrix_copy,
                                 .zero = petsc_matrix_zero,
                                 .scale = petsc_matrix_scale,
                                 .diag_scale = petsc_matrix_diag_scale,
                                 .add_identity = petsc_matrix_add_identity,
                                 .add_diagonal = petsc_matrix_add_diagonal,
                                 .set_diagonal = petsc_matrix_set_diagonal,
                                 .matvec = petsc_matrix_matvec,
                                 .set_values = petsc_matrix_set_values,
                                 .add_values = petsc_matrix_add_values,
                                 .get_values = petsc_matrix_get_values,
                                 .assemble = petsc_matrix_assemble,
                                 .dtor = petsc_matrix_dtor};
  return krylov_matrix_new(A, vtable, A->comm, (int)N_local, N_global);
}

static int petsc_matrix_block_size(void* context, index_t row)
{
  petsc_matrix_t* mat = context;
  if (mat->block_sizes != NULL)
    return mat->block_sizes[row];
  else
    return mat->block_size;
}

static void petsc_matrix_manipulate_fixed_blocks(void* context, index_t num_blocks,
                                                 index_t* block_rows, index_t* block_columns, 
                                                 real_t* block_values,
                                                 void (*manipulate)(void*, index_t, index_t*, index_t*, index_t*, real_t*),
                                                 bool copy_out)
{
  petsc_matrix_t* mat = context;
  int bs = mat->block_size;

  // We simply treat the blocks one at a time.
  for (index_t i = 0; i < num_blocks; ++i)
  {
    // Assemble the rows/columns.
    index_t block_row = block_rows[i];
    index_t block_column = block_columns[i];
    index_t num_rows = bs;
    index_t rows[bs], num_columns[bs], columns[bs*bs];
    {
      int l = 0;
      for (int j = 0; j < bs; ++j)
      {
        rows[j] = bs * block_row + j;
        num_columns[j] = bs;
        for (int k = 0; k < bs; ++k, ++l)
          columns[l] = bs * block_column + k;
      }
    }

    // Copy in the values if we are inserting/adding.
    real_t values[bs*bs];
    if (!copy_out)
    {
      int l = 0;
      for (int j = 0; j < bs; ++j)
      {
        rows[j] = bs * block_row + j;
        num_columns[j] = bs;
        for (int k = 0; k < bs; ++k, ++l)
          values[l] = block_values[k*bs+j];
      }
    }

    // Manipulate the values.
    manipulate(context, num_rows, num_columns, rows, columns, values);

    // Copy out the values if we are reading.
    if (copy_out)
    {
      int l = 0;
      for (int j = 0; j < bs; ++j)
      {
        rows[j] = bs * block_row + j;
        num_columns[j] = bs;
        for (int k = 0; k < bs; ++k, ++l)
          block_values[k*bs+j] = values[l];
      }
    }
  }
}

static void petsc_matrix_set_fixed_blocks(void* context, index_t num_blocks,
                                          index_t* block_rows, index_t* block_columns, 
                                          real_t* block_values)
{
  petsc_matrix_manipulate_fixed_blocks(context, num_blocks, 
                                       block_rows, block_columns, block_values,
                                       petsc_matrix_set_values, false);
}

static void petsc_matrix_add_fixed_blocks(void* context, index_t num_blocks,
                                          index_t* block_rows, index_t* block_columns, 
                                          real_t* block_values)
{
  petsc_matrix_manipulate_fixed_blocks(context, num_blocks, 
                                       block_rows, block_columns, block_values,
                                       petsc_matrix_add_values, false);
}

static void petsc_matrix_get_fixed_blocks(void* context, index_t num_blocks,
                                          index_t* block_rows, index_t* block_columns, 
                                          real_t* block_values)
{
  petsc_matrix_manipulate_fixed_blocks(context, num_blocks, 
                                       block_rows, block_columns, block_values,
                                       petsc_matrix_get_values, true);
}

static krylov_matrix_t* petsc_factory_block_matrix(void* context,
                                                   matrix_sparsity_t* sparsity,
                                                   int block_size)
{
  petsc_matrix_t* A = polymec_malloc(sizeof(petsc_matrix_t));
  A->factory = context;
  A->comm = matrix_sparsity_comm(sparsity);
  A->block_size = block_size;
  A->block_sizes = NULL;

  // PETSc deals in sparsities in terms of rows/columns, not blocks rows/columns.
  // So we have to create a block sparsity that contains the entire non-zero
  // structure of the matrix.
  matrix_sparsity_t* block_sp = matrix_sparsity_with_block_size(sparsity, block_size);
  index_t N_local = matrix_sparsity_num_local_rows(block_sp);
  index_t N_global = matrix_sparsity_num_global_rows(block_sp);

  int nprocs, rank;
  MPI_Comm_size(A->comm, &nprocs);
  MPI_Comm_rank(A->comm, &rank);
  A->factory->methods.MatCreate(A->comm, &A->A);
  if (nprocs == 1)
    A->factory->methods.MatSetType(A->A, MATSEQBAIJ);
  else
    A->factory->methods.MatSetType(A->A, MATMPIBAIJ);

  A->factory->methods.MatSetSizes(A->A, N_local, N_local, N_global, N_global);
  A->factory->methods.MatSetBlockSize(A->A, block_size);

  index_t* row_dist = matrix_sparsity_row_distribution(sparsity);
  PetscInt start = row_dist[rank], end = row_dist[rank+1];

  // Perform preallocation.
  // Note that nonzeros in preallocation are counted in blocks, not nonzeros.
  if (nprocs == 1)
  {
    PetscInt nnz[N_local];
    int pos = 0, r = 0;
    index_t row;
    while (matrix_sparsity_next_row(sparsity, &pos, &row))
    {
      nnz[r] = (PetscInt)(matrix_sparsity_num_columns(sparsity, row));
      ++r;
    }
    A->factory->methods.MatSeqBAIJSetPreallocation(A->A, (PetscInt)block_size, PETSC_DEFAULT, nnz);
  }
  else
  {
    // d_nnz[r] is the number of nonzeros for the "diagonal" portion of row r, 
    // and o_nnz[r] is the number of nonzeros for the "off-diagonal" portion.
    PetscInt d_nnz[N_local], o_nnz[N_local];
    memset(d_nnz, 0, sizeof(PetscInt) * N_local);
    memset(o_nnz, 0, sizeof(PetscInt) * N_local);
    index_t row, r = 0;
    int rpos = 0;
    while (matrix_sparsity_next_row(sparsity, &rpos, &row))
    {
      int cpos = 0;
      index_t column;
      while (matrix_sparsity_next_column(sparsity, row, &cpos, &column))
      {
        if ((column < start) || (column >= end))
          o_nnz[r] += 1;
        else
          d_nnz[r] += 1;
      }
      ++r;
    }
    A->factory->methods.MatMPIBAIJSetPreallocation(A->A, (PetscInt)block_size, PETSC_DEFAULT, d_nnz, PETSC_DEFAULT, o_nnz);
  }

  // Set the "non-zero" values to zero initially. This constructs the specific non-zero structure.
  real_t B[block_size*block_size];
  memset(B, 0, sizeof(real_t) * block_size*block_size);
  index_t row;
  int rpos = 0;
  while (matrix_sparsity_next_row(sparsity, &rpos, &row))
  {
    int cpos = 0;
    index_t column;
    while (matrix_sparsity_next_column(sparsity, row, &cpos, &column))
      petsc_matrix_set_fixed_blocks(A, 1, &row, &column, B);
  }

  // Assemble the matrix.
  petsc_matrix_assemble(A);
  matrix_sparsity_free(block_sp);

  // Set up the virtual table.
  krylov_matrix_vtable vtable = {.clone = petsc_matrix_clone,
                                 .copy = petsc_matrix_copy,
                                 .zero = petsc_matrix_zero,
                                 .scale = petsc_matrix_scale,
                                 .diag_scale = petsc_matrix_diag_scale,
                                 .add_identity = petsc_matrix_add_identity,
                                 .add_diagonal = petsc_matrix_add_diagonal,
                                 .set_diagonal = petsc_matrix_set_diagonal,
                                 .matvec = petsc_matrix_matvec,
                                 .set_values = petsc_matrix_set_values,
                                 .add_values = petsc_matrix_add_values,
                                 .get_values = petsc_matrix_get_values,
                                 .block_size = petsc_matrix_block_size,
                                 .set_blocks = petsc_matrix_set_fixed_blocks,
                                 .add_blocks = petsc_matrix_add_fixed_blocks,
                                 .get_blocks = petsc_matrix_get_fixed_blocks,
                                 .assemble = petsc_matrix_assemble,
                                 .fprintf = petsc_matrix_fprintf,
                                 .dtor = petsc_matrix_dtor};
  return krylov_matrix_new(A, vtable, A->comm, (int)N_local, N_global);
}

static void petsc_matrix_manipulate_var_blocks(void* context, index_t num_blocks,
                                               index_t* block_rows, index_t* block_columns, 
                                               real_t* block_values,
                                               void (*manipulate)(void*, index_t, index_t*, index_t*, index_t*, real_t*),
                                               bool copy_out)
{
  petsc_matrix_t* mat = context;

  // We simply treat the blocks one at a time.
  for (index_t i = 0; i < num_blocks; ++i)
  {
    // Assemble the rows/columns.
    index_t block_row = block_rows[i];
    index_t block_column = block_columns[i];
    int bs = petsc_matrix_block_size(mat, block_row);
    index_t num_rows = bs;
    index_t rows[bs], num_columns[bs], columns[bs*bs];
    {
      int l = 0;
      for (int j = 0; j < bs; ++j)
      {
        rows[j] = bs * block_row + j;
        num_columns[j] = bs;
        for (int k = 0; k < bs; ++k, ++l)
          columns[l] = bs * block_column + k;
      }
    }

    // Copy in the values if we are inserting/adding.
    real_t values[bs*bs];
    if (!copy_out)
    {
      int l = 0;
      for (int j = 0; j < bs; ++j)
      {
        rows[j] = bs * block_row + j;
        num_columns[j] = bs;
        for (int k = 0; k < bs; ++k, ++l)
          values[l] = block_values[k*bs+j];
      }
    }

    // Manipulate the values.
    manipulate(context, num_rows, num_columns, rows, columns, values);

    // Copy out the values if we are reading.
    if (copy_out)
    {
      int l = 0;
      for (int j = 0; j < bs; ++j)
      {
        rows[j] = bs * block_row + j;
        num_columns[j] = bs;
        for (int k = 0; k < bs; ++k, ++l)
          block_values[k*bs+j] = values[l];
      }
    }
  }
}

static void petsc_matrix_set_var_blocks(void* context, index_t num_blocks,
                                          index_t* block_rows, index_t* block_columns, 
                                          real_t* block_values)
{
  petsc_matrix_manipulate_var_blocks(context, num_blocks, 
                                     block_rows, block_columns, block_values,
                                     petsc_matrix_set_values, false);
}

static void petsc_matrix_add_var_blocks(void* context, index_t num_blocks,
                                        index_t* block_rows, index_t* block_columns, 
                                        real_t* block_values)
{
  petsc_matrix_manipulate_var_blocks(context, num_blocks, 
                                     block_rows, block_columns, block_values,
                                     petsc_matrix_add_values, false);
}

static void petsc_matrix_get_var_blocks(void* context, index_t num_blocks,
                                        index_t* block_rows, index_t* block_columns, 
                                        real_t* block_values)
{
  petsc_matrix_manipulate_var_blocks(context, num_blocks, 
                                     block_rows, block_columns, block_values,
                                     petsc_matrix_get_values, true);
}

// We replicate this struct here from krylov_solver.c in order to override
// methods in the matrix vtable to fake block matrices.
struct krylov_matrix_t
{
  void* context;
  krylov_matrix_vtable vtable;
  MPI_Comm comm;
  int num_local_rows;
  int num_global_rows;
};

static krylov_matrix_t* petsc_factory_var_block_matrix(void* context,
                                                       matrix_sparsity_t* sparsity,
                                                       int* block_sizes)
{
  // Create a block matrix sparsity pattern with the given block sizes, 
  // and make a regular matrix from that.
  matrix_sparsity_t* block_sp = matrix_sparsity_with_block_sizes(sparsity, block_sizes);
  krylov_matrix_t* mat = petsc_factory_matrix(context, block_sp);
  matrix_sparsity_free(block_sp);

  // Set the block size and override the block matrix methods.
  petsc_matrix_t* A = mat->context;
  index_t N_local = matrix_sparsity_num_local_rows(sparsity);
  A->block_sizes = polymec_malloc(sizeof(int) * N_local);
  memcpy(A->block_sizes, block_sizes, sizeof(int) * N_local);
  mat->vtable.block_size = petsc_matrix_block_size;
  mat->vtable.set_blocks = petsc_matrix_set_var_blocks;
  mat->vtable.add_blocks = petsc_matrix_add_var_blocks;
  mat->vtable.get_blocks = petsc_matrix_get_var_blocks;
  return mat;
}

static void* petsc_vector_clone(void* context)
{
  petsc_vector_t* v = context;
  petsc_vector_t* clone = polymec_malloc(sizeof(petsc_vector_t));
  clone->factory = v->factory;
  clone->comm = v->comm;
  PetscErrorCode err = v->factory->methods.VecDuplicate(v->v, &(clone->v));
  if (err != 0)
    polymec_error("petsc_vector_clone: Error in VecDuplicate!");
  err = v->factory->methods.VecCopy(v->v, clone->v);
  if (err != 0)
    polymec_error("petsc_vector_clone: Error in VecCopy!");
  return clone;
}

static void petsc_vector_copy(void* context, void* copy)
{
  petsc_vector_t* v = context;
  petsc_vector_t* v1 = copy;
  PetscErrorCode err = v->factory->methods.VecCopy(v->v, v1->v);
  if (err != 0)
    polymec_error("petsc_vector_clone: Error in VecCopy!");
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

static void petsc_vector_diag_scale(void* context, void* D)
{
  petsc_vector_t* v = context;
  petsc_vector_t* d = D;
  PetscInt size;
  v->factory->methods.VecGetLocalSize(v->v, &size);
  real_t* vi;
  v->factory->methods.VecGetArray(v->v, &vi);
  real_t* Di;
  d->factory->methods.VecGetArray(d->v, &Di);
  for (PetscInt i = 0; i < size; ++i)
    vi[i] *= Di[i];
  v->factory->methods.VecRestoreArray(v->v, &vi);
  d->factory->methods.VecRestoreArray(v->v, &Di);
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

static void petsc_vector_get_values(void* context, index_t num_values,
                                    index_t* indices, real_t* values)
{
  petsc_vector_t* v = context;
  v->factory->methods.VecGetValues(v->v, num_values, indices, values);
}

static void petsc_vector_copy_in(void* context, real_t* local_values)
{
  petsc_vector_t* v = context;
  PetscInt size;
  v->factory->methods.VecGetLocalSize(v->v, &size);
  real_t* array;
  v->factory->methods.VecGetArray(v->v, &array);
  memcpy(array, local_values, sizeof(real_t) * size);
  v->factory->methods.VecRestoreArray(v->v, &array);
}

static void petsc_vector_copy_out(void* context, real_t* local_values)
{
  petsc_vector_t* v = context;
  PetscInt size;
  v->factory->methods.VecGetLocalSize(v->v, &size);
  real_t* array;
  v->factory->methods.VecGetArray(v->v, &array);
  memcpy(local_values, array, sizeof(real_t) * size);
  v->factory->methods.VecRestoreArray(v->v, &array);
}

static real_t petsc_vector_dot(void* context, void* W)
{
  petsc_vector_t* v = context;
  petsc_vector_t* w = W;
  real_t dot;
  v->factory->methods.VecDot(v->v, w->v, &dot);
  return dot;
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

static real_t petsc_vector_w2_norm(void* context, void* W)
{
  petsc_vector_t* v = context;
  petsc_vector_t* w = W;

  // Accumulate the local part of the norm.
  PetscInt n;
  v->factory->methods.VecGetLocalSize(v->v, &n);
  real_t* v_values;
  v->factory->methods.VecGetArray(v->v, &v_values);
  real_t* w_values;
  w->factory->methods.VecGetArray(w->v, &w_values);
  real_t local_norm = 0.0;
  for (PetscInt i = 0; i < n; ++i)
  {
    real_t wi = w_values[i];
    real_t vi = v_values[i];
    local_norm += (real_t)(wi*wi*vi*vi);
  }
  v->factory->methods.VecRestoreArray(v->v, &v_values);
  w->factory->methods.VecRestoreArray(w->v, &w_values);

  // Now mash together all the parallel portions.
  real_t global_norm;
  MPI_Allreduce(&local_norm, &global_norm, 1, MPI_REAL_T, MPI_SUM, v->comm);
  return sqrt(global_norm);
}

static real_t petsc_vector_wrms_norm(void* context, void* W)
{
  petsc_vector_t* v = context;
  petsc_vector_t* w = W;

  // Accumulate the local part of the norm.
  PetscInt n;
  v->factory->methods.VecGetLocalSize(v->v, &n);
  real_t* v_values;
  v->factory->methods.VecGetArray(v->v, &v_values);
  real_t* w_values;
  w->factory->methods.VecGetArray(w->v, &w_values);
  real_t local_norm = 0.0;
  for (PetscInt i = 0; i < n; ++i)
  {
    real_t wi = w_values[i];
    real_t vi = v_values[i];
    local_norm += (real_t)(wi*wi*vi*vi);
  }
  v->factory->methods.VecRestoreArray(v->v, &v_values);
  w->factory->methods.VecRestoreArray(w->v, &w_values);

  // Now mash together all the parallel portions.
  real_t local_data[2] = {local_norm, (real_t)n};
  real_t global_data[2];
  MPI_Allreduce(local_data, global_data, 2, MPI_REAL_T, MPI_SUM, v->comm);
  return sqrt(global_data[0]/global_data[1]);
}

static void petsc_vector_fprintf(void* context, FILE* stream)
{
  petsc_vector_t* v = context;
  PetscViewer viewer;
  v->factory->methods.PetscViewerASCIIOpenWithFILE(v->comm, stream, &viewer);
  v->factory->methods.VecView(v->v, viewer);
  v->factory->methods.PetscViewerDestroy(&viewer);
}

static void petsc_vector_dtor(void* context)
{
  petsc_vector_t* v = context;
  v->factory->methods.VecDestroy(&(v->v));
  v->factory = NULL;
  polymec_free(v);
}

static krylov_vector_t* petsc_factory_vector(void* context,
                                             MPI_Comm comm,
                                             index_t* row_dist)
{
  petsc_vector_t* v = polymec_malloc(sizeof(petsc_vector_t));
  v->factory = context;
  v->comm = comm;

  int rank, nprocs;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &nprocs);

  index_t N_local = row_dist[rank+1] - row_dist[rank];
  index_t N_global = row_dist[nprocs];

  if (nprocs == 1)
    v->factory->methods.VecCreateSeq(comm, N_global, &v->v);
  else
    v->factory->methods.VecCreateMPI(comm, N_local, N_global, &v->v);

  petsc_vector_assemble(v);
  petsc_vector_zero(v);

  // Set up the virtual table.
  krylov_vector_vtable vtable = {.clone = petsc_vector_clone,
                                 .copy = petsc_vector_copy,
                                 .zero = petsc_vector_zero,
                                 .set_value = petsc_vector_set_value,
                                 .scale = petsc_vector_scale,
                                 .diag_scale = petsc_vector_diag_scale,
                                 .set_values = petsc_vector_set_values,
                                 .add_values = petsc_vector_add_values,
                                 .get_values = petsc_vector_get_values,
                                 .copy_in = petsc_vector_copy_in,
                                 .copy_out = petsc_vector_copy_out,
                                 .assemble = petsc_vector_assemble,
                                 .dot = petsc_vector_dot,
                                 .norm = petsc_vector_norm,
                                 .w2_norm = petsc_vector_w2_norm,
                                 .wrms_norm = petsc_vector_wrms_norm,
                                 .fprintf = petsc_vector_fprintf,
                                 .dtor = petsc_vector_dtor};
  return krylov_vector_new(v, vtable, (int)N_local, N_global);
}

static void petsc_factory_dtor(void* context)
{
  petsc_factory_t* factory = context;
  MPI_Barrier(MPI_COMM_WORLD);
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

// This handles PETSc errors as Polymec errors.
static PetscErrorCode polymec_petsc_error_handler(MPI_Comm comm,
                                                  int line,
                                                  const char* func,
                                                  const char* file,
                                                  PetscErrorCode n,
                                                  int p,
                                                  const char* mess,
                                                  void* ctx)
{
  polymec_error("PETSC error detected in %s (%s:%d): %s", func, file, line, mess);
}

#endif

krylov_factory_t* petsc_krylov_factory(const char* petsc_dir,
                                       const char* petsc_arch)
{
#if POLYMEC_HAVE_SHARED_LIBS
#if POLYMEC_HAVE_MPI
  ASSERT(((petsc_dir == NULL) && (petsc_arch == NULL)) ||
         ((petsc_dir != NULL) && (petsc_arch != NULL)));

  petsc_factory_t* factory = polymec_malloc(sizeof(petsc_factory_t));

  // Try to find PETSc.
  char petsc_path[FILENAME_MAX+1];
  char *my_petsc_dir = NULL, *my_petsc_arch = NULL;
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
    int pos = 0;
    size_t length;
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

  // Get the symbols.
#define FETCH_PETSC_SYMBOL(symbol_name) \
  FETCH_SYMBOL(petsc, #symbol_name, factory->methods.symbol_name, failure);

  FETCH_PETSC_SYMBOL(PetscInitializeNoArguments);
  FETCH_PETSC_SYMBOL(PetscInitialized);
  FETCH_PETSC_SYMBOL(PetscFinalize);
  FETCH_PETSC_SYMBOL(PetscPushErrorHandler);
  FETCH_PETSC_SYMBOL(PETSC_VIEWER_STDOUT_);
  FETCH_PETSC_SYMBOL(PetscViewerASCIIOpenWithFILE);
  FETCH_PETSC_SYMBOL(PetscViewerDestroy);

  FETCH_PETSC_SYMBOL(KSPCreate);
  FETCH_PETSC_SYMBOL(KSPSetType);
  FETCH_PETSC_SYMBOL(KSPSetFromOptions);
  FETCH_PETSC_SYMBOL(KSPGetPC);
  FETCH_PETSC_SYMBOL(KSPSetPC);
  FETCH_PETSC_SYMBOL(KSPSetTolerances);
  FETCH_PETSC_SYMBOL(KSPGetTolerances);
  FETCH_PETSC_SYMBOL(KSPSetDiagonalScale);
  FETCH_PETSC_SYMBOL(KSPSetDiagonalScaleFix);
  FETCH_PETSC_SYMBOL(KSPSetNormType);
  FETCH_PETSC_SYMBOL(KSPSetOperators);
  FETCH_PETSC_SYMBOL(KSPSetUp);
  FETCH_PETSC_SYMBOL(KSPSolve);
  FETCH_PETSC_SYMBOL(KSPGetConvergedReason);
  FETCH_PETSC_SYMBOL(KSPGetIterationNumber);
  FETCH_PETSC_SYMBOL(KSPGetResidualNorm);
  FETCH_PETSC_SYMBOL(KSPDestroy);
  FETCH_PETSC_SYMBOL(KSPGMRESSetRestart);

  FETCH_PETSC_SYMBOL(PCCreate);
  FETCH_PETSC_SYMBOL(PCSetType);
  FETCH_PETSC_SYMBOL(PCSetFromOptions);

  FETCH_PETSC_SYMBOL(MatCreate);
  FETCH_PETSC_SYMBOL(MatConvert);
  FETCH_PETSC_SYMBOL(MatCopy);
  FETCH_PETSC_SYMBOL(MatSetType);
  FETCH_PETSC_SYMBOL(MatSetSizes);
  FETCH_PETSC_SYMBOL(MatGetOwnershipRange);
  FETCH_PETSC_SYMBOL(MatSeqAIJSetPreallocation);
  FETCH_PETSC_SYMBOL(MatMPIAIJSetPreallocation);
  FETCH_PETSC_SYMBOL(MatSeqBAIJSetPreallocation);
  FETCH_PETSC_SYMBOL(MatMPIBAIJSetPreallocation);
  FETCH_PETSC_SYMBOL(MatSetBlockSize);
  FETCH_PETSC_SYMBOL(MatDestroy);
  FETCH_PETSC_SYMBOL(MatScale);
  FETCH_PETSC_SYMBOL(MatDiagonalScale);
  FETCH_PETSC_SYMBOL(MatShift);
  FETCH_PETSC_SYMBOL(MatDiagonalSet);
  FETCH_PETSC_SYMBOL(MatZeroEntries);
  FETCH_PETSC_SYMBOL(MatMult);
  FETCH_PETSC_SYMBOL(MatMultTranspose);
  FETCH_PETSC_SYMBOL(MatGetSize);
  FETCH_PETSC_SYMBOL(MatGetLocalSize);
  FETCH_PETSC_SYMBOL(MatSetValues);
  FETCH_PETSC_SYMBOL(MatGetValues);
  FETCH_PETSC_SYMBOL(MatAssemblyBegin);
  FETCH_PETSC_SYMBOL(MatAssemblyEnd);
  FETCH_PETSC_SYMBOL(MatView);
  
  FETCH_PETSC_SYMBOL(VecCreateSeq);
  FETCH_PETSC_SYMBOL(VecCreateMPI);
  FETCH_PETSC_SYMBOL(VecDuplicate);
  FETCH_PETSC_SYMBOL(VecCopy);
  FETCH_PETSC_SYMBOL(VecSetUp);
  FETCH_PETSC_SYMBOL(VecDestroy);
  FETCH_PETSC_SYMBOL(VecGetSize);
  FETCH_PETSC_SYMBOL(VecGetLocalSize);
  FETCH_PETSC_SYMBOL(VecZeroEntries);
  FETCH_PETSC_SYMBOL(VecScale);
  FETCH_PETSC_SYMBOL(VecSet);
  FETCH_PETSC_SYMBOL(VecSetValues);
  FETCH_PETSC_SYMBOL(VecGetValues);
  FETCH_PETSC_SYMBOL(VecGetArray);
  FETCH_PETSC_SYMBOL(VecRestoreArray);
  FETCH_PETSC_SYMBOL(VecAssemblyBegin);
  FETCH_PETSC_SYMBOL(VecAssemblyEnd);
  FETCH_PETSC_SYMBOL(VecDot);
  FETCH_PETSC_SYMBOL(VecNorm);
  FETCH_PETSC_SYMBOL(VecView);
#undef FETCH_PETSC_SYMBOL

  log_debug("petsc_krylov_factory: Got PETSc symbols.");

  // Initialize PETSc if needed.
  factory->finalize_petsc = false;
  PetscBool initialized;
  factory->methods.PetscInitialized(&initialized);
  if (initialized == PETSC_FALSE)
  {
    log_debug("petsc_krylov_factory: Initializing PETSc.");
    factory->methods.PetscInitializeNoArguments();

    // Since we are the first to use this PETSc library, we must finalize it.
    factory->finalize_petsc = true;
  }

  // Hook into polymec's error handler.
  factory->methods.PetscPushErrorHandler(polymec_petsc_error_handler, NULL);

  // Stash the library.
  factory->petsc = petsc; 

  // Construct the factory.
  krylov_factory_vtable vtable = {.pcg_solver = petsc_factory_pcg_solver,
                                  .gmres_solver = petsc_factory_gmres_solver,
                                  .bicgstab_solver = petsc_factory_bicgstab_solver,
                                  .special_solver = petsc_factory_special_solver,
                                  .preconditioner = petsc_factory_pc,
                                  .matrix = petsc_factory_matrix,
                                  .block_matrix = petsc_factory_block_matrix,
                                  .var_block_matrix = petsc_factory_var_block_matrix,
                                  .vector = petsc_factory_vector,
                                  .dtor = petsc_factory_dtor};
  char name[1024];
  snprintf(name, 1023, "PETSc v%s", version);
  return krylov_factory_new(name, factory, vtable);

failure:
  dlclose(petsc);
  polymec_free(factory);
  return NULL;
#else
  log_urgent("petsc_krylov_factory: Polymec must be configured with MPI to use PETSc.");
  return NULL;
#endif
#else
  log_urgent("petsc_krylov_factory: Polymec must be configured with MPI and shared library support to use PETSc.");
  return NULL;
#endif
}

