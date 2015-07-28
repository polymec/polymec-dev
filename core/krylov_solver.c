// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <float.h>
#include <dlfcn.h>
#include "core/krylov_solver.h"
#include "core/timer.h"

//------------------------------------------------------------------------
//                Krylov data types and virtual tables
//------------------------------------------------------------------------

typedef struct
{
  void (*set_tolerances)(void* context, real_t rel_tol, real_t abs_tol, real_t div_tol);
  void (*set_max_iterations)(void* context, int max_iters);
  void (*set_operator)(void* context, void* op);
  bool (*solve)(void* context, void* x, void* b, real_t* res_norm, int* num_iters);
  void (*dtor)(void* context);
} krylov_solver_vtable;

struct krylov_solver_t
{
  char* name;
  void* context;
  krylov_solver_vtable vtable;
};

typedef struct
{
  void* (*clone)(void* context);
  void (*zero)(void* context);
  void (*scale)(void* context, real_t scale_factor);
  void (*add_identity)(void* context, real_t scale_factor);
  void (*add_diagonal)(void* context, void* D);
  void (*set_diagonal)(void* context, void* D);
  void (*set_values)(void* context, int num_rows, int* num_columns, index_t* rows, index_t* columns, real_t* values);
  void (*add_values)(void* context, int num_rows, int* num_columns, index_t* rows, index_t* columns, real_t* values);
  void (*start_assembly)(void* context);
  void (*finish_assembly)(void* context);
  void (*get_values)(void* context, int num_rows, int* num_columns, index_t* rows, index_t* columns, real_t* values);
  void (*dtor)(void* context);
} krylov_matrix_vtable;

struct krylov_matrix_t
{
  void* context;
  krylov_matrix_vtable vtable;
  int num_local_rows, num_global_rows;
  int num_local_columns, num_global_columns;
};

typedef struct
{
  void* (*clone)(void* context);
  void (*zero)(void* context);
  void (*set_value)(void* context, real_t value);
  void (*scale)(void* context, real_t scale_factor);
  void (*set_values)(void* context, int num_values, index_t* indices, real_t* values);
  void (*add_values)(void* context, int num_values, index_t* indices, real_t* values);
  void (*start_assembly)(void* context);
  void (*finish_assembly)(void* context);
  void (*get_values)(void* context, int num_values, index_t* indices, real_t* values);
  real_t (*norm)(void* context, int p);
  void (*dtor)(void* context);
} krylov_vector_vtable;

struct krylov_vector_t
{
  void* context;
  krylov_vector_vtable vtable;
  int local_size, global_size;
};

typedef struct 
{
  krylov_solver_t* (*solver)(void* context, MPI_Comm comm, string_string_unordered_map_t* options);
  krylov_matrix_t* (*matrix)(void* context, adj_graph_t* sparsity);
  krylov_matrix_t* (*block_matrix)(void* context, adj_graph_t* sparsity, int block_size);
  krylov_vector_t* (*vector)(void* context, MPI_Comm comm, int N);
  void (*dtor)(void* context);
} krylov_factory_vtable;

struct krylov_factory_t
{
  // Parallel stuff.
  int rank, nprocs;
  MPI_Comm comm;

  char* name;
  void* context;
  krylov_factory_vtable vtable;
};


//------------------------------------------------------------------------
//                          Krylov solver
//------------------------------------------------------------------------

static krylov_solver_t* krylov_solver_new(const char* name,
                                          void* context,
                                          krylov_solver_vtable vtable)
{
  krylov_solver_t* solver = polymec_malloc(sizeof(krylov_solver_t));
  solver->name = string_dup(name);
  solver->context = context;
  solver->vtable = vtable;
  return solver;
}

char* krylov_solver_name(krylov_solver_t* solver)
{
  return solver->name;
}

void krylov_solver_free(krylov_solver_t* solver)
{
  if ((solver->context != NULL) && (solver->vtable.dtor != NULL))
    solver->vtable.dtor(solver->context);
  string_free(solver->name);
  polymec_free(solver);
}

void* krylov_solver_impl(krylov_solver_t* solver)
{
  return solver->context;
}

void krylov_solver_set_tolerances(krylov_solver_t* solver, 
                                  real_t relative_tolerance,
                                  real_t absolute_tolerance,
                                  real_t divergence_tolerance)
{
  ASSERT(relative_tolerance > 0.0);
  ASSERT(absolute_tolerance > 0.0);
  ASSERT(divergence_tolerance > 0.0);
  solver->vtable.set_tolerances(solver->context, 
                                relative_tolerance,
                                absolute_tolerance,
                                divergence_tolerance);
}

void krylov_solver_set_max_iterations(krylov_solver_t* solver, 
                                      int max_iterations)
{
  ASSERT(max_iterations > 0);
  solver->vtable.set_max_iterations(solver->context, max_iterations);
}

void krylov_solver_set_operator(krylov_solver_t* solver, 
                                krylov_matrix_t* op)
{
  solver->vtable.set_operator(solver->context, op->context);
}

bool krylov_solver_solve(krylov_solver_t* solver, 
                         krylov_vector_t* x, 
                         krylov_vector_t* b, 
                         real_t* residual_norm, 
                         int* num_iterations)
{
  return solver->vtable.solve(solver->context, x->context, b->context,
                              residual_norm, num_iterations);
}

//------------------------------------------------------------------------
//                          Krylov matrix
//------------------------------------------------------------------------

static krylov_matrix_t* krylov_matrix_new(void* context,
                                          krylov_matrix_vtable vtable)
{
  ASSERT(vtable.zero != NULL);
  ASSERT(vtable.add_diagonal != NULL);
  ASSERT(vtable.set_diagonal != NULL);
  ASSERT(vtable.set_values != NULL);
  ASSERT(vtable.add_values != NULL);
  ASSERT(vtable.get_values != NULL);
  krylov_matrix_t* A = polymec_malloc(sizeof(krylov_matrix_t));
  A->context = context;
  A->vtable = vtable;
  return A;
}

void krylov_matrix_free(krylov_matrix_t* A)
{
  if ((A->context != NULL) && (A->vtable.dtor != NULL))
    A->vtable.dtor(A->context);
  polymec_free(A);
}

krylov_matrix_t* krylov_matrix_clone(krylov_matrix_t* A)
{
  krylov_matrix_t* B = polymec_malloc(sizeof(krylov_matrix_t));
  B->context = A->vtable.clone(A->context);
  B->vtable = A->vtable;
  B->num_local_rows = A->num_local_rows;
  B->num_local_columns = A->num_local_columns;
  B->num_global_rows = A->num_global_rows;
  B->num_global_columns = A->num_global_columns;
  return B;
}

void* krylov_matrix_impl(krylov_matrix_t* A)
{
  return A->context;
}

int krylov_matrix_num_local_rows(krylov_matrix_t* A)
{
  return A->num_local_rows;
}

int krylov_matrix_num_local_columns(krylov_matrix_t* A)
{
  return A->num_local_columns;
}

int krylov_matrix_num_global_rows(krylov_matrix_t* A)
{
  return A->num_global_rows;
}

int krylov_matrix_num_global_columns(krylov_matrix_t* A)
{
  return A->num_global_columns;
}

void krylov_matrix_zero(krylov_matrix_t* A)
{
  A->vtable.zero(A->context);
}

void krylov_matrix_add_identity(krylov_matrix_t* A,
                                real_t scale_factor)
{
  A->vtable.add_identity(A->context, scale_factor);
}

void krylov_matrix_scale(krylov_matrix_t* A,
                         real_t scale_factor)
{
  A->vtable.scale(A->context, scale_factor);
}

void krylov_matrix_add_diagonal(krylov_matrix_t* A,
                                krylov_matrix_t* D)
{
  A->vtable.add_diagonal(A->context, D->context);
}

void krylov_matrix_set_diagonal(krylov_matrix_t* A,
                                krylov_matrix_t* D)
{
  A->vtable.set_diagonal(A->context, D->context);
}

void krylov_matrix_set_values(krylov_matrix_t* A,
                              int num_rows,
                              int* num_columns,
                              index_t* rows, index_t* columns,
                              real_t* values)
{
  A->vtable.set_values(A->context, num_rows, num_columns, rows, columns, values);
}
                              
void krylov_matrix_add_values(krylov_matrix_t* A,
                              int num_rows,
                              int* num_columns,
                              index_t* rows, index_t* columns,
                              real_t* values)
{
  A->vtable.add_values(A->context, num_rows, num_columns, rows, columns, values);
}
                              
void krylov_matrix_assemble(krylov_matrix_t* A)
{
  krylov_matrix_start_assembly(A);
  krylov_matrix_finish_assembly(A);
}

void krylov_matrix_start_assembly(krylov_matrix_t* A)
{
  if (A->vtable.start_assembly != NULL);
    A->vtable.start_assembly(A->context);
}

void krylov_matrix_finish_assembly(krylov_matrix_t* A)
{
  if (A->vtable.finish_assembly != NULL);
    A->vtable.finish_assembly(A->context);
}

void krylov_matrix_get_values(krylov_matrix_t* A,
                              int num_rows,
                              int* num_columns,
                              index_t* rows, index_t* columns,
                              real_t* values)
{
  A->vtable.get_values(A->context, num_rows, num_columns, rows, columns, values);
}

//------------------------------------------------------------------------
//                          Krylov vector
//------------------------------------------------------------------------

static krylov_vector_t* krylov_vector_new(void* context,
                                          krylov_vector_vtable vtable)
{
  ASSERT(vtable.zero != NULL);
  ASSERT(vtable.set_value != NULL);
  ASSERT(vtable.scale != NULL);
  ASSERT(vtable.set_values != NULL);
  ASSERT(vtable.add_values != NULL);
  ASSERT(vtable.get_values != NULL);
  krylov_vector_t* v = polymec_malloc(sizeof(krylov_vector_t));
  v->context = context;
  v->vtable = vtable;
  return v;
}

void krylov_vector_free(krylov_vector_t* v)
{
  if ((v->context != NULL) && (v->vtable.dtor != NULL))
    v->vtable.dtor(v->context);
  polymec_free(v);
}

krylov_vector_t* krylov_vector_clone(krylov_vector_t* v)
{
  krylov_vector_t* u = polymec_malloc(sizeof(krylov_vector_t));
  u->context = v->vtable.clone(v->context);
  u->vtable = v->vtable;
  u->local_size = v->local_size;
  u->global_size = v->global_size;
  return u;
}

void* krylov_vector_impl(krylov_vector_t* v)
{
  return v->context;
}

int krylov_vector_local_size(krylov_vector_t* v)
{
  return v->local_size;
}

int krylov_vector_global_size(krylov_vector_t* v)
{
  return v->global_size;
}

void krylov_vector_zero(krylov_vector_t* v)
{
  v->vtable.zero(v->context);
}

void krylov_vector_set_value(krylov_vector_t* v,
                             real_t value)
{
  v->vtable.set_value(v->context, value);
}

void krylov_vector_scale(krylov_vector_t* v,
                         real_t scale_factor)
{
  v->vtable.scale(v->context, scale_factor);
}

void krylov_vector_set_values(krylov_vector_t* v,
                              int num_values,
                              index_t* indices,
                              real_t* values)
{
  v->vtable.set_values(v->context, num_values, indices, values);
}
                              
void krylov_vector_add_values(krylov_vector_t* v,
                              int num_values,
                              index_t* indices,
                              real_t* values)
{
  v->vtable.add_values(v->context, num_values, indices, values);
}

void krylov_vector_assemble(krylov_vector_t* v)
{
  krylov_vector_start_assembly(v);
  krylov_vector_finish_assembly(v);
}

void krylov_vector_start_assembly(krylov_vector_t* v)
{
  if (v->vtable.start_assembly != NULL);
    v->vtable.start_assembly(v->context);
}

void krylov_vector_finish_assembly(krylov_vector_t* v)
{
  if (v->vtable.finish_assembly != NULL);
    v->vtable.finish_assembly(v->context);
}

void krylov_vector_get_values(krylov_vector_t* v,
                              int num_values,
                              index_t* indices,
                              real_t* values)
{
  v->vtable.get_values(v->context, num_values, indices, values);
}

real_t krylov_vector_norm(krylov_vector_t* v, int p)
{
  ASSERT((p == 0) || (p == 1) || (p == 2));
  return v->vtable.norm(v->context, p);
}

//------------------------------------------------------------------------
//                  Factories for creating Krylov solvers
//------------------------------------------------------------------------

static krylov_factory_t* krylov_factory_new(const char* name,
                                            void* context,
                                            krylov_factory_vtable vtable)
{
  ASSERT(vtable.solver != NULL);
  ASSERT(vtable.matrix != NULL);
  ASSERT(vtable.block_matrix != NULL);
  ASSERT(vtable.vector != NULL);
  krylov_factory_t* factory = polymec_malloc(sizeof(krylov_factory_t));
  factory->name = string_dup(name);
  factory->context = context;
  factory->vtable = vtable;
  return factory;
}

void krylov_factory_free(krylov_factory_t* factory)
{
  if ((factory->context != NULL) && (factory->vtable.dtor != NULL))
    factory->vtable.dtor(factory->context);
  string_free(factory->name);
  polymec_free(factory);
}

char* krylov_factory_name(krylov_factory_t* factory)
{
  return factory->name;
}

krylov_matrix_t* krylov_factory_matrix(krylov_factory_t* factory, 
                                       adj_graph_t* sparsity)
{
  return factory->vtable.matrix(factory->context, sparsity);
}

krylov_matrix_t* krylov_factory_block_matrix(krylov_factory_t* factory, 
                                             adj_graph_t* sparsity,
                                             int block_size)
{
  ASSERT(block_size > 0);
  return factory->vtable.block_matrix(factory->context, sparsity, block_size);
}

krylov_vector_t* krylov_factory_vector(krylov_factory_t* factory,
                                       MPI_Comm comm,
                                       int N)
{
  ASSERT(N > 0);
  return factory->vtable.vector(factory->context, comm, N);
}

krylov_solver_t* krylov_factory_solver(krylov_factory_t* factory,
                                       MPI_Comm comm,
                                       string_string_unordered_map_t* options)
{
  return factory->vtable.solver(factory->context, comm, options);
}

//------------------------------------------------------------------------
//                          PETSc implementation
//------------------------------------------------------------------------

// Here's a table of function pointers for the PETSc library.
typedef real_t PetscScalar;
typedef real_t PetscReal;
typedef index_t PetscInt;
typedef int PetscErrorCode;
typedef void* KSP;
typedef const char* KSPType;
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
#define PETSC_DECIDE -1
#define PETSC_DETERMINE PETSC_DECIDE
typedef struct
{
  PetscErrorCode (*KSPCreate)(MPI_Comm,KSP *);
  PetscErrorCode (*KSPSetType)(KSP,KSPType);
  PetscErrorCode (*KSPSetUp)(KSP);
  PetscErrorCode (*KSPSetTolerances)(KSP,PetscReal,PetscReal,PetscReal,PetscInt);
  PetscErrorCode (*KSPGetTolerances)(KSP,PetscReal*,PetscReal*,PetscReal*,PetscInt*);
  PetscErrorCode (*KSPSetOperators)(KSP,Mat,Mat);
  PetscErrorCode (*KSPSolve)(KSP,Vec,Vec);
  PetscErrorCode (*KSPGetConvergedReason)(KSP,int*);
  PetscErrorCode (*KSPGetIterationNumber)(KSP,PetscInt*);
  PetscErrorCode (*KSPGetResidualNorm)(KSP,PetscReal*);
  PetscErrorCode (*KSPDestroy)(KSP*);
  PetscErrorCode (*MatCreate)(MPI_Comm,Mat*);
  PetscErrorCode (*MatConvert)(Mat, MatType, MatReuse, Mat*);
  PetscErrorCode (*MatSetType)(Mat, MatType);
  PetscErrorCode (*MatSetSizes)(Mat, PetscInt m, PetscInt n, PetscInt M, PetscInt N);
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
  petsc_methods_table methods;
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
  polymec_free(solver);
}

static krylov_solver_t* petsc_factory_solver(void* context,
                                             MPI_Comm comm,
                                             string_string_unordered_map_t* options)
{
  petsc_solver_t* solver = polymec_malloc(sizeof(petsc_solver_t));
  solver->factory = context;
  solver->factory->methods.KSPCreate(comm, &solver->ksp);
  char** type_p = string_string_unordered_map_get(options, "type");
  if (type_p == NULL)
    solver->factory->methods.KSPSetType(solver->ksp, "gmres");
  else
    solver->factory->methods.KSPSetType(solver->ksp, *type_p);
  // Set up the virtual table.
  krylov_solver_vtable vtable = {.set_tolerances = petsc_solver_set_tolerances,
                                 .set_max_iterations = petsc_solver_set_max_iterations,
                                 .set_operator = petsc_solver_set_operator,
                                 .solve = petsc_solver_solve,
                                 .dtor = petsc_solver_dtor};
  return krylov_solver_new("KSP", solver, vtable);
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
  A->factory->methods.MatDiagonalSet(A->A, diag, INSERT_VALUES);
}

static void petsc_matrix_enter_values(void* context, int num_rows,
                                      int* num_columns, index_t* rows, index_t* columns,
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

static void petsc_matrix_set_values(void* context, int num_rows,
                                    int* num_columns, index_t* rows, index_t* columns,
                                    real_t* values)
{
  petsc_matrix_enter_values(context, num_rows, num_columns, rows, columns, values, INSERT_VALUES);
}

static void petsc_matrix_add_values(void* context, int num_rows,
                                    int* num_columns, index_t* rows, index_t* columns,
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

static void petsc_matrix_get_values(void* context, int num_rows,
                                    int* num_columns, index_t* rows, index_t* columns,
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
  // FIXME: Do preallocation here!

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
  return krylov_matrix_new(A, vtable);
}

static krylov_matrix_t* petsc_factory_block_matrix(void* context,
                                                   adj_graph_t* sparsity,
                                                   int block_size)
{
  petsc_matrix_t* A = polymec_malloc(sizeof(petsc_matrix_t));
  A->factory = context;
  // FIXME
  // Set up the virtual table.
  krylov_matrix_vtable vtable = {.clone = petsc_matrix_clone,
                                 .zero = petsc_matrix_zero,
                                 .scale = petsc_matrix_scale,
                                 .add_identity = petsc_matrix_add_identity,
                                 .set_diagonal = petsc_matrix_set_diagonal,
                                 .set_values = petsc_matrix_set_values,
                                 .add_values = petsc_matrix_add_values,
                                 .start_assembly = petsc_matrix_start_assembly,
                                 .finish_assembly = petsc_matrix_finish_assembly,
                                 .get_values = petsc_matrix_get_values,
                                 .dtor = petsc_matrix_dtor};
  return krylov_matrix_new(A, vtable);
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

static void petsc_vector_set_values(void* context, int num_values,
                                    index_t* indices, real_t* values)
{
  petsc_vector_t* v = context;
  v->factory->methods.VecSetValues(v->v, num_values, indices, values, INSERT_VALUES);
}

static void petsc_vector_add_values(void* context, int num_values,
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

static void petsc_vector_get_values(void* context, int num_values,
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
  polymec_free(v);
}

static krylov_vector_t* petsc_factory_vector(void* context,
                                             MPI_Comm comm,
                                             int N)
{
  petsc_vector_t* v = polymec_malloc(sizeof(petsc_vector_t));
  v->factory = context;
  int nprocs;
  MPI_Comm_size(comm, &nprocs);
  if (nprocs == 1)
    v->factory->methods.VecCreateSeq(comm, N, &v->v);
  else
    v->factory->methods.VecCreateMPI(comm, N, PETSC_DETERMINE, &v->v);
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
  return krylov_vector_new(v, vtable);
}

krylov_factory_t* petsc_krylov_factory(const char* petsc_dir,
                                       const char* petsc_arch)
{
  ASSERT(((petsc_dir == NULL) && (petsc_arch == NULL)) ||
          (petsc_dir != NULL) && (petsc_arch != NULL));

  char name[128];
  petsc_factory_t* factory = polymec_malloc(sizeof(petsc_factory_t));

  // Try to find PETSc.
  char petsc_path[FILENAME_MAX+1];
#ifdef APPLE
  const char* dl_suffix = "dylib";
#else
  const char* dl_suffix = "so";
#endif
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
    snprintf(petsc_path, FILENAME_MAX, "%s/%s/lib/libpetsc.%s", petsc_dir, petsc_arch, dl_suffix);
  else
    snprintf(petsc_path, FILENAME_MAX, "%s/%s/lib/libpetsc.%s", my_petsc_dir, my_petsc_arch, dl_suffix);

  // Try to open libPETSc and mine it for symbols.
  void* petsc = dlopen(petsc_path, RTLD_NOW);
  if (petsc == NULL)
  {
    char* msg = dlerror();
    polymec_error("petsc_krylov_factory: %s.", msg);
  }

  // Non-standard C (but POSIX compliant!). Thanks, dlsym.
  *((void**)&(factory->methods.KSPCreate)) = dlsym(petsc, "KSPCreate");
  *((void**)&(factory->methods.KSPSetType)) = dlsym(petsc, "KSPSetType");
  *((void**)&(factory->methods.KSPSetUp)) = dlsym(petsc, "KSPSetUp");
  *((void**)&(factory->methods.KSPSetTolerances)) = dlsym(petsc, "KSPSetTolerances");
  *((void**)&(factory->methods.KSPGetTolerances)) = dlsym(petsc, "KSPGetTolerances");
  *((void**)&(factory->methods.KSPSetOperators)) = dlsym(petsc, "KSPSetOperators");
  *((void**)&(factory->methods.KSPSolve)) = dlsym(petsc, "KSPSolve");
  *((void**)&(factory->methods.KSPGetIterationNumber)) = dlsym(petsc, "KSPGetIterationNumber");
  *((void**)&(factory->methods.KSPGetResidualNorm)) = dlsym(petsc, "KSPGetResidualNorm");
  *((void**)&(factory->methods.KSPDestroy)) = dlsym(petsc, "KSPDestroy");

  *((void**)&(factory->methods.MatCreate)) = dlsym(petsc, "MatCreate");
  *((void**)&(factory->methods.MatConvert)) = dlsym(petsc, "MatConvert");
  *((void**)&(factory->methods.MatSetType)) = dlsym(petsc, "MatSetType");
  *((void**)&(factory->methods.MatSetSizes)) = dlsym(petsc, "MatSetSizes");
  *((void**)&(factory->methods.MatSetUp)) = dlsym(petsc, "MatSetUp");
  *((void**)&(factory->methods.MatDestroy)) = dlsym(petsc, "MatDestroy");
  *((void**)&(factory->methods.MatScale)) = dlsym(petsc, "MatScale");
  *((void**)&(factory->methods.MatShift)) = dlsym(petsc, "MatShift");
  *((void**)&(factory->methods.MatDiagonalSet)) = dlsym(petsc, "MatDiagonalSet");
  *((void**)&(factory->methods.MatZeroEntries)) = dlsym(petsc, "MatZeroEntries");
  *((void**)&(factory->methods.MatGetSize)) = dlsym(petsc, "MatGetSize");
  *((void**)&(factory->methods.MatGetLocalSize)) = dlsym(petsc, "MatGetLocalSize");
  *((void**)&(factory->methods.MatSetValues)) = dlsym(petsc, "MatSetValues");
  *((void**)&(factory->methods.MatGetValues)) = dlsym(petsc, "MatGetValues");
  *((void**)&(factory->methods.MatAssemblyBegin)) = dlsym(petsc, "MatAssemblyBegin");
  *((void**)&(factory->methods.MatAssemblyEnd)) = dlsym(petsc, "MatAssemblyEnd");
  
  *((void**)&(factory->methods.VecCreateSeq)) = dlsym(petsc, "VecCreateSeq");
  *((void**)&(factory->methods.VecCreateMPI)) = dlsym(petsc, "VecCreateMPI");
  *((void**)&(factory->methods.VecDuplicate)) = dlsym(petsc, "VecDuplicate");
  *((void**)&(factory->methods.VecCopy)) = dlsym(petsc, "VecCopy");
  *((void**)&(factory->methods.VecSetUp)) = dlsym(petsc, "VecSetUp");
  *((void**)&(factory->methods.VecDestroy)) = dlsym(petsc, "VecDestroy");
  *((void**)&(factory->methods.VecZeroEntries)) = dlsym(petsc, "VecZeroEntries");
  *((void**)&(factory->methods.VecScale)) = dlsym(petsc, "VecScale");
  *((void**)&(factory->methods.VecSet)) = dlsym(petsc, "VecSet");
  *((void**)&(factory->methods.VecSetValues)) = dlsym(petsc, "VecSetValues");
  *((void**)&(factory->methods.VecGetValues)) = dlsym(petsc, "VecGetValues");
  *((void**)&(factory->methods.VecAssemblyBegin)) = dlsym(petsc, "VecAssemblyBegin");
  *((void**)&(factory->methods.VecAssemblyEnd)) = dlsym(petsc, "VecAssemblyEnd");
  *((void**)&(factory->methods.VecNorm)) = dlsym(petsc, "VecNorm");

  // Finish up and construct the factory.
  dlclose(petsc);
  krylov_factory_vtable vtable = {.solver = petsc_factory_solver,
                                  .matrix = petsc_factory_matrix,
                                  .block_matrix = petsc_factory_block_matrix,
                                  .vector = petsc_factory_vector};
  return krylov_factory_new(name, factory, vtable);
}

//------------------------------------------------------------------------
//                          HYPRE implementation
//------------------------------------------------------------------------

typedef struct
{
  int i;
} hypre_methods_table;

typedef struct
{
  hypre_methods_table methods;
} hypre_factory_t;

krylov_factory_t* hypre_krylov_factory(const char* library_path)
{
  POLYMEC_NOT_IMPLEMENTED;
//  char name[128];
//  hypre_factory_t* factory = polymec_malloc(sizeof(hypre_factory_t));
//  krylov_factory_vtable vtable;
//  return krylov_factory_new(name, factory, vtable);
}

