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
  void (*set_values)(void* context, int num_rows, int* num_columns, int* rows, int* columns, real_t* values);
  void (*add_values)(void* context, int num_rows, int* num_columns, int* rows, int* columns, real_t* values);
  void (*start_assembly)(void* context);
  void (*finish_assembly)(void* context);
  void (*get_values)(void* context, int num_rows, int* num_columns, int* rows, int* columns, real_t* values);
  void (*dtor)(void* context);
} krylov_matrix_vtable;

struct krylov_matrix_t
{
  void* context;
  krylov_matrix_vtable vtable;
  int num_local_rows;
  int num_global_rows;
};

typedef struct
{
  void* (*clone)(void* context);
  void (*zero)(void* context);
  void (*set_value)(void* context, real_t value);
  void (*scale)(void* context, real_t scale_factor);
  void (*set_values)(void* context, int num_values, int* indices, real_t* values);
  void (*add_values)(void* context, int num_values, int* indices, real_t* values);
  void (*start_assembly)(void* context);
  void (*finish_assembly)(void* context);
  void (*get_values)(void* context, int num_values, int* indices, real_t* values);
  real_t (*norm)(void* context, int p);
  void (*dtor)(void* context);
} krylov_vector_vtable;

struct krylov_vector_t
{
  void* context;
  krylov_vector_vtable vtable;
  int local_size;
  int global_size;
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

void krylov_solver_free(krylov_solver_t* solver)
{
  if ((solver->context != NULL) && (solver->vtable.dtor != NULL))
    solver->vtable.dtor(solver->context);
  string_free(solver->name);
  polymec_free(solver);
}

char* krylov_solver_name(krylov_solver_t* solver)
{
  return solver->name;
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
                                          krylov_matrix_vtable vtable,
                                          int num_local_rows,
                                          int num_global_rows)
{
  ASSERT(vtable.zero != NULL);
  ASSERT(vtable.add_diagonal != NULL);
  ASSERT(vtable.set_diagonal != NULL);
  ASSERT(vtable.set_values != NULL);
  ASSERT(vtable.add_values != NULL);
  ASSERT(vtable.get_values != NULL);
  ASSERT(num_local_rows > 0);
  ASSERT(num_global_rows > 0);
  krylov_matrix_t* A = polymec_malloc(sizeof(krylov_matrix_t));
  A->context = context;
  A->vtable = vtable;
  A->num_local_rows = num_local_rows;
  A->num_global_rows = num_global_rows;
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
  B->num_global_rows = A->num_global_rows;
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

int krylov_matrix_num_global_rows(krylov_matrix_t* A)
{
  return A->num_global_rows;
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
                              int* rows, int* columns,
                              real_t* values)
{
  A->vtable.set_values(A->context, num_rows, num_columns, rows, columns, values);
}
                              
void krylov_matrix_add_values(krylov_matrix_t* A,
                              int num_rows,
                              int* num_columns,
                              int* rows, int* columns,
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
                              int* rows, int* columns,
                              real_t* values)
{
  A->vtable.get_values(A->context, num_rows, num_columns, rows, columns, values);
}

//------------------------------------------------------------------------
//                          Krylov vector
//------------------------------------------------------------------------

static krylov_vector_t* krylov_vector_new(void* context,
                                          krylov_vector_vtable vtable,
                                          int local_size,
                                          int global_size)
{
  ASSERT(vtable.zero != NULL);
  ASSERT(vtable.set_value != NULL);
  ASSERT(vtable.scale != NULL);
  ASSERT(vtable.set_values != NULL);
  ASSERT(vtable.add_values != NULL);
  ASSERT(vtable.get_values != NULL);
  ASSERT(local_size > 0);
  ASSERT(global_size > 0);
  krylov_vector_t* v = polymec_malloc(sizeof(krylov_vector_t));
  v->context = context;
  v->vtable = vtable;
  v->local_size = local_size;
  v->global_size = global_size;
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
                              int* indices,
                              real_t* values)
{
  v->vtable.set_values(v->context, num_values, indices, values);
}
                              
void krylov_vector_add_values(krylov_vector_t* v,
                              int num_values,
                              int* indices,
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
                              int* indices,
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

//------------------------------------------------------------------------
//                          PETSc implementation
//------------------------------------------------------------------------

// Here's a table of function pointers for the PETSc library.
typedef real_t PetscScalar;
typedef real_t PetscReal;
typedef int PetscMPIInt;
typedef int PetscInt; // 32-bit indices for PETSc.
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
  PetscErrorCode (*PCSetFromOptions)(PC);
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

static krylov_solver_t* petsc_factory_solver(void* context,
                                             MPI_Comm comm,
                                             string_string_unordered_map_t* options)
{
  petsc_solver_t* solver = polymec_malloc(sizeof(petsc_solver_t));
  solver->factory = context;
  solver->factory->methods.KSPCreate(comm, &solver->ksp);

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
                                      int* num_columns, int* rows, int* columns,
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
                                    int* num_columns, int* rows, int* columns,
                                    real_t* values)
{
  petsc_matrix_enter_values(context, num_rows, num_columns, rows, columns, values, INSERT_VALUES);
}

static void petsc_matrix_add_values(void* context, int num_rows,
                                    int* num_columns, int* rows, int* columns,
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
                                    int* num_columns, int* rows, int* columns,
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

static void petsc_vector_set_values(void* context, int num_values,
                                    int* indices, real_t* values)
{
  petsc_vector_t* v = context;
  v->factory->methods.VecSetValues(v->v, num_values, indices, values, INSERT_VALUES);
}

static void petsc_vector_add_values(void* context, int num_values,
                                    int* indices, real_t* values)
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
                                    int* indices, real_t* values)
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
  PetscInt N_global;
  v->factory->methods.VecGetSize(v->v, &N_global);
printf("Global size: %zd\n", N_global);
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
  return krylov_vector_new(v, vtable, N, N_global);
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

krylov_factory_t* petsc_krylov_factory(const char* petsc_dir,
                                       const char* petsc_arch)
{
  ASSERT(((petsc_dir == NULL) && (petsc_arch == NULL)) ||
          (petsc_dir != NULL) && (petsc_arch != NULL));

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
      string_copy_from_raw(line, length+1, text);
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

  // Get the other symbols.
  FETCH_SYMBOL(petsc, "PetscInitialize", factory->methods.PetscInitialize, failure);
  FETCH_SYMBOL(petsc, "PetscInitialized", factory->methods.PetscInitialized, failure);
  FETCH_SYMBOL(petsc, "PetscFinalize", factory->methods.PetscFinalize, failure);

  FETCH_SYMBOL(petsc, "KSPCreate", factory->methods.KSPCreate, failure);
  FETCH_SYMBOL(petsc, "KSPSetType", factory->methods.KSPSetType, failure);
  FETCH_SYMBOL(petsc, "KSPSetFromOptions", factory->methods.KSPSetFromOptions, failure);
  FETCH_SYMBOL(petsc, "KSPSetUp", factory->methods.KSPSetUp, failure);
  FETCH_SYMBOL(petsc, "KSPGetPC", factory->methods.KSPGetPC, failure);
  FETCH_SYMBOL(petsc, "PCSetFromOptions", factory->methods.PCSetFromOptions, failure);
  FETCH_SYMBOL(petsc, "KSPSetTolerances", factory->methods.KSPSetTolerances, failure);
  FETCH_SYMBOL(petsc, "KSPGetTolerances", factory->methods.KSPGetTolerances, failure);
  FETCH_SYMBOL(petsc, "KSPSetOperators", factory->methods.KSPSetOperators, failure);
  FETCH_SYMBOL(petsc, "KSPSolve", factory->methods.KSPSolve, failure);
  FETCH_SYMBOL(petsc, "KSPGetConvergedReason", factory->methods.KSPGetConvergedReason, failure);
  FETCH_SYMBOL(petsc, "KSPGetIterationNumber", factory->methods.KSPGetIterationNumber, failure);
  FETCH_SYMBOL(petsc, "KSPGetResidualNorm", factory->methods.KSPGetResidualNorm, failure);
  FETCH_SYMBOL(petsc, "KSPDestroy", factory->methods.KSPDestroy, failure);

  FETCH_SYMBOL(petsc, "MatCreate", factory->methods.MatCreate, failure);
  FETCH_SYMBOL(petsc, "MatConvert", factory->methods.MatConvert, failure);
  FETCH_SYMBOL(petsc, "MatSetType", factory->methods.MatSetType, failure);
  FETCH_SYMBOL(petsc, "MatSetSizes", factory->methods.MatSetSizes, failure);
  FETCH_SYMBOL(petsc, "MatSeqAIJSetPreallocation", factory->methods.MatSeqAIJSetPreallocation, failure);
  FETCH_SYMBOL(petsc, "MatMPIAIJSetPreallocation", factory->methods.MatMPIAIJSetPreallocation, failure);
  FETCH_SYMBOL(petsc, "MatSeqBAIJSetPreallocation", factory->methods.MatSeqBAIJSetPreallocation, failure);
  FETCH_SYMBOL(petsc, "MatMPIBAIJSetPreallocation", factory->methods.MatMPIBAIJSetPreallocation, failure);
  FETCH_SYMBOL(petsc, "MatSetBlockSize", factory->methods.MatSetBlockSize, failure);
  FETCH_SYMBOL(petsc, "MatSetUp", factory->methods.MatSetUp, failure);
  FETCH_SYMBOL(petsc, "MatDestroy", factory->methods.MatDestroy, failure);
  FETCH_SYMBOL(petsc, "MatScale", factory->methods.MatScale, failure);
  FETCH_SYMBOL(petsc, "MatShift", factory->methods.MatShift, failure);
  FETCH_SYMBOL(petsc, "MatDiagonalSet", factory->methods.MatDiagonalSet, failure);
  FETCH_SYMBOL(petsc, "MatZeroEntries", factory->methods.MatZeroEntries, failure);
  FETCH_SYMBOL(petsc, "MatGetSize", factory->methods.MatGetSize, failure);
  FETCH_SYMBOL(petsc, "MatGetLocalSize", factory->methods.MatGetLocalSize, failure);
  FETCH_SYMBOL(petsc, "MatSetValues", factory->methods.MatSetValues, failure);
  FETCH_SYMBOL(petsc, "MatGetValues", factory->methods.MatGetValues, failure);
  FETCH_SYMBOL(petsc, "MatAssemblyBegin", factory->methods.MatAssemblyBegin, failure);
  FETCH_SYMBOL(petsc, "MatAssemblyEnd", factory->methods.MatAssemblyEnd, failure);
  
  FETCH_SYMBOL(petsc, "VecCreateSeq", factory->methods.VecCreateSeq, failure);
  FETCH_SYMBOL(petsc, "VecCreateMPI", factory->methods.VecCreateMPI, failure);
  FETCH_SYMBOL(petsc, "VecDuplicate", factory->methods.VecDuplicate, failure);
  FETCH_SYMBOL(petsc, "VecCopy", factory->methods.VecCopy, failure);
  FETCH_SYMBOL(petsc, "VecSetUp", factory->methods.VecSetUp, failure);
  FETCH_SYMBOL(petsc, "VecDestroy", factory->methods.VecDestroy, failure);
  FETCH_SYMBOL(petsc, "VecGetSize", factory->methods.VecGetSize, failure);
  FETCH_SYMBOL(petsc, "VecZeroEntries", factory->methods.VecZeroEntries, failure);
  FETCH_SYMBOL(petsc, "VecScale", factory->methods.VecScale, failure);
  FETCH_SYMBOL(petsc, "VecSet", factory->methods.VecSet, failure);
  FETCH_SYMBOL(petsc, "VecSetValues", factory->methods.VecSetValues, failure);
  FETCH_SYMBOL(petsc, "VecGetValues", factory->methods.VecGetValues, failure);
  FETCH_SYMBOL(petsc, "VecAssemblyBegin", factory->methods.VecAssemblyBegin, failure);
  FETCH_SYMBOL(petsc, "VecAssemblyEnd", factory->methods.VecAssemblyEnd, failure);
  FETCH_SYMBOL(petsc, "VecNorm", factory->methods.VecNorm, failure);

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
  krylov_factory_vtable vtable = {.solver = petsc_factory_solver,
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

