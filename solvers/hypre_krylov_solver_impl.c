// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// This file gets included from elsewhere, and requires the following 
// things to be defined:
// HypreFactory - A string containing the name of the factory constructor.
// HYPRE_Int - The integer type used for indexing.
// HYPRE_Real - The type used to store real floating point values.
// HYPRE_Complex - The type used to store complex floating point values.

#include <float.h>
#include <dlfcn.h>
#include "core/polymec.h"
#include "core/timer.h"
#include "core/text_buffer.h"
#include "core/slist.h"
#include "core/array.h"
#include "solvers/krylov_solver.h"

// HYPRE types and definitions.
typedef real_t HYPRE_Real;
typedef complex_t HYPRE_Complex;
#define HYPRE_PARCSR  5555
typedef void* HYPRE_Solver;
typedef void* HYPRE_Matrix;
typedef void* HYPRE_Vector;
typedef void* HYPRE_IJMatrix;
typedef void* HYPRE_IJVector;
typedef void* HYPRE_ParCSRMatrix;
typedef void* HYPRE_ParVector;
typedef HYPRE_Int (*HYPRE_PtrToSolverFcn)(HYPRE_Solver, HYPRE_Matrix, HYPRE_Vector, HYPRE_Vector);
typedef HYPRE_Int (*HYPRE_PtrToParSolverFcn)(HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector);
#if POLYMEC_HAVE_MPI
typedef MPI_Comm HYPRE_MPI_Comm;
#else
typedef HYPRE_Int HYPRE_MPI_Comm;
#endif

// Types of supported solver methods.
typedef enum
{
  HYPRE_PCG,
  HYPRE_GMRES,
  HYPRE_BICGSTAB,
  HYPRE_BOOMERANG,
} hypre_solver_method_t;

// Here's a table of function pointers for the HYPRE library.
typedef struct
{
  HYPRE_Int (*HYPRE_GetError)(void);
  void      (*HYPRE_DescribeError)(HYPRE_Int, char*);
  HYPRE_Int (*HYPRE_ClearAllErrors)(void);

  HYPRE_Int (*HYPRE_ParCSRPCGCreate)(HYPRE_MPI_Comm, HYPRE_Solver*);
  HYPRE_Int (*HYPRE_ParCSRPCGDestroy)(HYPRE_Solver);
  HYPRE_Int (*HYPRE_ParCSRPCGSetup)(HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector);
  HYPRE_Int (*HYPRE_ParCSRPCGSolve)(HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector);
  HYPRE_Int (*HYPRE_ParCSRPCGSetTol)(HYPRE_Solver, HYPRE_Real);
  HYPRE_Int (*HYPRE_ParCSRPCGSetAbsoluteTol)(HYPRE_Solver, HYPRE_Real);
  HYPRE_Int (*HYPRE_ParCSRPCGSetMaxIter)(HYPRE_Solver, HYPRE_Int);
  HYPRE_Int (*HYPRE_ParCSRPCGSetStopCrit)(HYPRE_Solver, HYPRE_Int);
  HYPRE_Int (*HYPRE_PCGSetTwoNorm)(HYPRE_Solver, HYPRE_Int);
  HYPRE_Int (*HYPRE_ParCSRPCGSetPrecond)(HYPRE_Solver, HYPRE_PtrToParSolverFcn, HYPRE_PtrToParSolverFcn, HYPRE_Solver);
  HYPRE_Int (*HYPRE_ParCSRPCGGetPrecond)(HYPRE_Solver, HYPRE_Solver*);
  HYPRE_Int (*HYPRE_ParCSRPCGGetNumIterations)(HYPRE_Solver, HYPRE_Int*);
  HYPRE_Int (*HYPRE_ParCSRPCGGetFinalRelativeResidualNorm)(HYPRE_Solver, HYPRE_Real*);
  HYPRE_Int (*HYPRE_ParCSRPCGSetPrintLevel)(HYPRE_Solver, HYPRE_Int);

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
  HYPRE_Int (*HYPRE_ParCSRGMRESSetPrintLevel)(HYPRE_Solver, HYPRE_Int);

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
  HYPRE_Int (*HYPRE_ParCSRBiCGSTABSetPrintLevel)(HYPRE_Solver, HYPRE_Int);

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
  HYPRE_Int (*HYPRE_EuclidSolve)(HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector);
  HYPRE_Int (*HYPRE_EuclidSetLevel)(HYPRE_Solver, HYPRE_Int);
  HYPRE_Int (*HYPRE_EuclidSetBJ)(HYPRE_Solver, HYPRE_Int);
  HYPRE_Int (*HYPRE_EuclidSetSparseA)(HYPRE_Solver, HYPRE_Real);
  HYPRE_Int (*HYPRE_EuclidSetRowScale)(HYPRE_Solver, HYPRE_Int);
  HYPRE_Int (*HYPRE_EuclidSetILUT)(HYPRE_Solver, HYPRE_Real);

  HYPRE_Int (*HYPRE_ParaSailsCreate)(HYPRE_MPI_Comm, HYPRE_Solver*);
  HYPRE_Int (*HYPRE_ParaSailsDestroy)(HYPRE_Solver); 
  HYPRE_Int (*HYPRE_ParaSailsSetup)(HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector);
  HYPRE_Int (*HYPRE_ParaSailsSolve)(HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector);
  HYPRE_Int (*HYPRE_ParaSailsSetParams)(HYPRE_Solver, HYPRE_Real, HYPRE_Int);
  HYPRE_Int (*HYPRE_ParaSailsSetFilter)(HYPRE_Solver, HYPRE_Real);
  HYPRE_Int (*HYPRE_ParaSailsSetSym)(HYPRE_Solver, HYPRE_Int);
  HYPRE_Int (*HYPRE_ParaSailsSetLoadbal)(HYPRE_Solver, HYPRE_Real);
  HYPRE_Int (*HYPRE_ParaSailsSetReuse)(HYPRE_Solver, HYPRE_Int);

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
  HYPRE_Int (*HYPRE_IJMatrixGetRowCounts)(HYPRE_IJMatrix, HYPRE_Int, HYPRE_Int*, HYPRE_Int*);
  HYPRE_Int (*HYPRE_IJMatrixPrint)(HYPRE_IJMatrix, const char*);
  HYPRE_Int (*HYPRE_IJMatrixSetPrintLevel)(HYPRE_IJMatrix, HYPRE_Int);

  HYPRE_Int (*HYPRE_IJVectorCreate)(HYPRE_MPI_Comm, HYPRE_Int, HYPRE_Int, HYPRE_IJVector*);
  HYPRE_Int (*HYPRE_IJVectorDestroy)(HYPRE_IJVector);
  HYPRE_Int (*HYPRE_IJVectorInitialize)(HYPRE_IJVector);
  HYPRE_Int (*HYPRE_IJVectorSetValues)(HYPRE_IJVector, HYPRE_Int, const HYPRE_Int*, const HYPRE_Real*);
  HYPRE_Int (*HYPRE_IJVectorAddToValues)(HYPRE_IJVector, HYPRE_Int, const HYPRE_Int*, const HYPRE_Real*);
  HYPRE_Int (*HYPRE_IJVectorAssemble)(HYPRE_IJVector);
  HYPRE_Int (*HYPRE_IJVectorGetValues)(HYPRE_IJVector, HYPRE_Int, HYPRE_Int*, HYPRE_Real*);
  HYPRE_Int (*HYPRE_IJVectorSetObjectType)(HYPRE_IJVector, HYPRE_Int);
  HYPRE_Int (*HYPRE_IJVectorGetObject)(HYPRE_IJVector, void**);
  HYPRE_Int (*HYPRE_IJVectorPrint)(HYPRE_IJVector, const char*);

  // Boneheaded built-in preconditioner.
  HYPRE_Int (*HYPRE_ParCSRDiagScaleSetup)(HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector);
  HYPRE_Int (*HYPRE_ParCSRDiagScale)(HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector);

  // Internals needed for scaling and zeroing matrices.
  HYPRE_Int (*HYPRE_ParCSRMatrixGetRow)(void*, HYPRE_Int, HYPRE_Int*, HYPRE_Int**, HYPRE_Real**);
  HYPRE_Int (*HYPRE_ParCSRMatrixRestoreRow)(void*, HYPRE_Int, HYPRE_Int*, HYPRE_Int**, HYPRE_Real**);

  // Matrix-vector product.
  HYPRE_Int (*HYPRE_ParCSRMatrixMatvec)(HYPRE_Real, HYPRE_ParCSRMatrix*, HYPRE_ParVector, HYPRE_Real, HYPRE_ParVector*);
  HYPRE_Int (*HYPRE_ParCSRMatrixMatvecT)(HYPRE_Real, HYPRE_ParCSRMatrix*, HYPRE_ParVector, HYPRE_Real, HYPRE_ParVector*);
} hypre_methods_table;

// "HYPRE extension" methods -- from an auxiliary library.
typedef struct
{
  HYPRE_Int (*HYPRE_IJMatrixDiagScale)(HYPRE_IJMatrix, HYPRE_IJVector, HYPRE_IJVector);
} hypre_ext_methods_table;

typedef struct
{
  ptr_array_t* hypre_libs;
  string_array_t* hypre_lib_names;
  hypre_methods_table methods;

  void* hypre_ext;
  hypre_ext_methods_table ext_methods;
} hypre_factory_t;

typedef struct
{
  hypre_factory_t* factory;
  hypre_solver_method_t type;
  HYPRE_Solver solver;
  HYPRE_ParCSRMatrix op;
} hypre_solver_t;

typedef struct
{
  HYPRE_Solver pc;
  HYPRE_Int (*setup)(HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector);
  HYPRE_Int (*solve)(HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector);
  HYPRE_Int (*dtor)(HYPRE_Solver);
} hypre_pc_t;

typedef struct
{
  MPI_Comm comm;
  HYPRE_Int ilow, ihigh;
  hypre_factory_t* factory;
  HYPRE_IJMatrix A;
  size_t block_size;
  size_t* block_sizes;
} hypre_matrix_t;

typedef struct
{
  MPI_Comm comm;
  HYPRE_Int ilow, ihigh;
  hypre_factory_t* factory;
  HYPRE_IJVector v;
} hypre_vector_t;

static bool solver_error_occurred(hypre_solver_t* s, 
                                  const char* func_name)
{
  HYPRE_Int error = s->factory->methods.HYPRE_GetError();
  if (error != 0)
  {
    char msg[1024];
    s->factory->methods.HYPRE_DescribeError(error, msg);
    s->factory->methods.HYPRE_ClearAllErrors();
    log_urgent("%s: %s", func_name, msg);
  }
  return (error != 0);
}

#define SOLVER_ERROR_OCCURRED(s) (solver_error_occurred(s, __func__))

static void hypre_solver_set_tolerances(void* context,
                                        real_t rel_tol,
                                        real_t abs_tol,
                                        real_t div_tol)

{
  hypre_solver_t* solver = context;

  // Dispatch.
  HYPRE_Int (*set_tol)(HYPRE_Solver, HYPRE_Real);
  HYPRE_Int (*set_abs_tol)(HYPRE_Solver, HYPRE_Real);
  switch(solver->type)
  {
    case HYPRE_PCG: 
      set_tol = solver->factory->methods.HYPRE_ParCSRPCGSetTol; 
      set_abs_tol = solver->factory->methods.HYPRE_ParCSRPCGSetAbsoluteTol; 
      break;
    case HYPRE_GMRES: 
      set_tol = solver->factory->methods.HYPRE_ParCSRGMRESSetTol; 
      set_abs_tol = solver->factory->methods.HYPRE_ParCSRGMRESSetAbsoluteTol; 
      break;
    case HYPRE_BICGSTAB:
      set_tol = solver->factory->methods.HYPRE_ParCSRBiCGSTABSetTol; 
      set_abs_tol = solver->factory->methods.HYPRE_ParCSRBiCGSTABSetAbsoluteTol; 
      break;
    case HYPRE_BOOMERANG:
      set_tol = solver->factory->methods.HYPRE_BoomerAMGSetTol; 
      set_abs_tol = NULL;
  }
  set_tol(solver->solver, (HYPRE_Real)rel_tol);
  if (set_abs_tol != NULL)
    set_abs_tol(solver->solver, (HYPRE_Real)abs_tol);
}

static void hypre_solver_set_max_iterations(void* context,
                                            int max_iters)
{
  hypre_solver_t* solver = context;

  // Dispatch.
  HYPRE_Int (*set_max_iters)(HYPRE_Solver, HYPRE_Int);
  switch(solver->type)
  {
    case HYPRE_PCG: 
      set_max_iters = solver->factory->methods.HYPRE_ParCSRPCGSetMaxIter; 
      break;
    case HYPRE_GMRES: 
      set_max_iters = solver->factory->methods.HYPRE_ParCSRGMRESSetMaxIter; 
      break;
    case HYPRE_BICGSTAB:
      set_max_iters = solver->factory->methods.HYPRE_ParCSRBiCGSTABSetMaxIter; 
      break;
    case HYPRE_BOOMERANG:
      set_max_iters = solver->factory->methods.HYPRE_BoomerAMGSetMaxIter; 
  }
  set_max_iters(solver->solver, (HYPRE_Int)max_iters);
}

static void hypre_solver_set_operator(void* context,
                                      void* op)
{
  hypre_solver_t* solver = context;
  HYPRE_IJMatrix A = ((hypre_matrix_t*)op)->A;
  solver->factory->methods.HYPRE_IJMatrixGetObject(A, &solver->op);

  // NOTE: We can't do HYPRE's full setup here, since we don't have 
  // NOTE: x or b.
}

static void hypre_solver_set_pc(void* context,
                                void* pc)
{
  hypre_solver_t* solver = context;
  hypre_pc_t* p = pc;

  // Dispatch.
  HYPRE_Int (*set_pc)(HYPRE_Solver, HYPRE_PtrToSolverFcn, HYPRE_PtrToSolverFcn, HYPRE_Solver);
  switch(solver->type)
  {
    case HYPRE_PCG: 
      set_pc = solver->factory->methods.HYPRE_ParCSRPCGSetPrecond; 
      break;
    case HYPRE_GMRES: 
      set_pc = solver->factory->methods.HYPRE_ParCSRGMRESSetPrecond; 
      break;
    case HYPRE_BICGSTAB:
      set_pc = solver->factory->methods.HYPRE_ParCSRBiCGSTABSetPrecond; 
      break;
    default:
      set_pc = NULL;
  }
  if (set_pc != NULL)
    set_pc(solver, (HYPRE_PtrToSolverFcn)p->solve, (HYPRE_PtrToSolverFcn)p->setup, p->pc);
}

static real_t hypre_vector_norm(void* context, int p);
static bool hypre_solver_solve(void* context,
                               void* b,
                               void* x,
                               real_t* res_norm,
                               int* num_iters)
{
  hypre_solver_t* solver = context;

  // Get our objects.
  HYPRE_IJVector B = ((hypre_vector_t*)b)->v;
  HYPRE_IJVector X = ((hypre_vector_t*)x)->v;
  HYPRE_ParVector par_B, par_X;
  solver->factory->methods.HYPRE_IJVectorGetObject(B, &par_B);
  solver->factory->methods.HYPRE_IJVectorGetObject(X, &par_X);

  // Record ||b||_2.
  real_t b_norm = hypre_vector_norm(b, 2);

  // Dispatch.
  HYPRE_Int (*setup)(HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector);
  HYPRE_Int (*solve)(HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector);
  HYPRE_Int (*get_num_iters)(HYPRE_Solver, HYPRE_Int*);
  HYPRE_Int (*get_norm)(HYPRE_Solver, HYPRE_Real*);
  switch(solver->type)
  {
    case HYPRE_PCG: 
      setup = solver->factory->methods.HYPRE_ParCSRPCGSetup; 
      solve = solver->factory->methods.HYPRE_ParCSRPCGSolve; 
      get_num_iters = solver->factory->methods.HYPRE_ParCSRPCGGetNumIterations;
      get_norm = solver->factory->methods.HYPRE_ParCSRPCGGetFinalRelativeResidualNorm;
      break;
    case HYPRE_GMRES: 
      setup = solver->factory->methods.HYPRE_ParCSRGMRESSetup; 
      solve = solver->factory->methods.HYPRE_ParCSRGMRESSolve; 
      get_num_iters = solver->factory->methods.HYPRE_ParCSRGMRESGetNumIterations;
      get_norm = solver->factory->methods.HYPRE_ParCSRGMRESGetFinalRelativeResidualNorm;
      break;
    case HYPRE_BICGSTAB:
      setup = solver->factory->methods.HYPRE_ParCSRBiCGSTABSetup; 
      solve = solver->factory->methods.HYPRE_ParCSRBiCGSTABSolve; 
      get_num_iters = solver->factory->methods.HYPRE_ParCSRBiCGSTABGetNumIterations; 
      get_norm = solver->factory->methods.HYPRE_ParCSRBiCGSTABGetFinalRelativeResidualNorm; 
      break;
    case HYPRE_BOOMERANG:
      setup = solver->factory->methods.HYPRE_BoomerAMGSetup; 
      solve = solver->factory->methods.HYPRE_BoomerAMGSolve; 
      get_num_iters = solver->factory->methods.HYPRE_BoomerAMGGetNumIterations; 
      get_norm = solver->factory->methods.HYPRE_BoomerAMGGetFinalRelativeResidualNorm; 
  }
  setup(solver->solver, solver->op, par_B, par_X);
  if (SOLVER_ERROR_OCCURRED(solver))
    return false;

  // HYPRE's Krylov methods can produce some pretty small numbers, which can be 
  // interpreted as denormalized garbage by polymec, so we tell polymec to chill.
  polymec_suspend_fpe();
  HYPRE_Int result = solve(solver->solver, solver->op, par_B, par_X);
  polymec_restore_fpe();

  int success = ((result == 0) && !SOLVER_ERROR_OCCURRED(solver));

  // Get the number of iterations.
  HYPRE_Int iters;
  get_num_iters(solver->solver, &iters);
  *num_iters = (int)iters;

  // Get the residual norm. HYPRE returns the relative norm, so we 
  // have to multiply by ||b||_2.
  HYPRE_Real norm;
  get_norm(solver->solver, &norm);
  *res_norm = norm * b_norm;

  solver->factory->methods.HYPRE_ClearAllErrors();
  return success;
}

static void hypre_solver_dtor(void* context)
{
  hypre_solver_t* solver = context;

  // Dispatch.
  HYPRE_Int (*destroy)(HYPRE_Solver);
  switch(solver->type)
  {
    case HYPRE_PCG: 
      destroy = solver->factory->methods.HYPRE_ParCSRPCGDestroy; 
      break;
    case HYPRE_GMRES: 
      destroy = solver->factory->methods.HYPRE_ParCSRGMRESDestroy; 
      break;
    case HYPRE_BICGSTAB:
      destroy = solver->factory->methods.HYPRE_ParCSRBiCGSTABDestroy; 
      break;
    case HYPRE_BOOMERANG:
      destroy = solver->factory->methods.HYPRE_BoomerAMGDestroy; 
  }
  destroy(solver->solver);
  solver->factory = NULL;
  polymec_free(solver);
}

static krylov_solver_t* hypre_factory_pcg_solver(void* context,
                                                 MPI_Comm comm)
     
{
  hypre_solver_t* solver = polymec_malloc(sizeof(hypre_solver_t));
  solver->factory = context;
  solver->type = HYPRE_PCG;
  HYPRE_MPI_Comm hypre_comm = comm;
  solver->factory->methods.HYPRE_ParCSRPCGCreate(hypre_comm, &solver->solver);
  if (SOLVER_ERROR_OCCURRED(solver))
  {
    polymec_free(solver);
    return NULL;
  }
  log_debug("hypre_factory_pcg_solver: Created PCG solver");

  // Set up diagonal scaling till someone tells us otherwise.
  solver->factory->methods.HYPRE_ParCSRPCGSetPrecond(solver->solver, 
                                                     (HYPRE_PtrToSolverFcn)solver->factory->methods.HYPRE_ParCSRDiagScale, 
                                                     (HYPRE_PtrToSolverFcn)solver->factory->methods.HYPRE_ParCSRDiagScaleSetup, 
                                                     solver->solver);
  if (SOLVER_ERROR_OCCURRED(solver))
  {
    polymec_free(solver);
    return NULL;
  }

  // We use the 2-norm for stopping criteria.
  solver->factory->methods.HYPRE_PCGSetTwoNorm(solver->solver, 1);

  // Set up the virtual table.
  krylov_solver_vtable vtable = {.set_tolerances = hypre_solver_set_tolerances,
                                 .set_max_iterations = hypre_solver_set_max_iterations,
                                 .set_operator = hypre_solver_set_operator,
                                 .set_preconditioner = hypre_solver_set_pc,
                                 .solve = hypre_solver_solve,
                                 .dtor = hypre_solver_dtor};
  return krylov_solver_new("HYPRE PCG", solver, vtable);
}

static krylov_solver_t* hypre_factory_gmres_solver(void* context,
                                                   MPI_Comm comm,
                                                   int krylov_dimension)
     
{
  hypre_solver_t* solver = polymec_malloc(sizeof(hypre_solver_t));
  solver->factory = context;
  solver->type = HYPRE_GMRES;
  HYPRE_MPI_Comm hypre_comm = comm;
  solver->factory->methods.HYPRE_ParCSRGMRESCreate(hypre_comm, &solver->solver);
  if (SOLVER_ERROR_OCCURRED(solver))
  {
    polymec_free(solver);
    return NULL;
  }

  solver->factory->methods.HYPRE_ParCSRGMRESSetKDim(solver->solver, (HYPRE_Int)krylov_dimension);
  if (SOLVER_ERROR_OCCURRED(solver))
  {
    polymec_free(solver);
    return NULL;
  }
  log_debug("hypre_factory_gmres_solver: Created solver with Krylov dim = %d", krylov_dimension);

  // Set up diagonal scaling till someone tells us otherwise.
  solver->factory->methods.HYPRE_ParCSRGMRESSetPrecond(solver->solver, 
                                                       (HYPRE_PtrToSolverFcn)solver->factory->methods.HYPRE_ParCSRDiagScale, 
                                                       (HYPRE_PtrToSolverFcn)solver->factory->methods.HYPRE_ParCSRDiagScaleSetup, 
                                                       solver->solver);
  if (SOLVER_ERROR_OCCURRED(solver))
  {
    polymec_free(solver);
    return NULL;
  }

  // Set up the virtual table.
  krylov_solver_vtable vtable = {.set_tolerances = hypre_solver_set_tolerances,
                                 .set_max_iterations = hypre_solver_set_max_iterations,
                                 .set_operator = hypre_solver_set_operator,
                                 .set_preconditioner = hypre_solver_set_pc,
                                 .solve = hypre_solver_solve,
                                 .dtor = hypre_solver_dtor};
  return krylov_solver_new("HYPRE GMRES", solver, vtable);
}

static krylov_solver_t* hypre_factory_bicgstab_solver(void* context,
                                                      MPI_Comm comm)
{
  hypre_solver_t* solver = polymec_malloc(sizeof(hypre_solver_t));
  solver->factory = context;
  solver->type = HYPRE_BICGSTAB;
  HYPRE_MPI_Comm hypre_comm = comm;
  solver->factory->methods.HYPRE_ParCSRBiCGSTABCreate(hypre_comm, &solver->solver);
  if (SOLVER_ERROR_OCCURRED(solver))
  {
    polymec_free(solver);
    return NULL;
  }

  // Set up diagonal scaling till someone tells us otherwise.
  solver->factory->methods.HYPRE_ParCSRBiCGSTABSetPrecond(solver->solver, 
                                                          (HYPRE_PtrToSolverFcn)solver->factory->methods.HYPRE_ParCSRDiagScale, 
                                                          (HYPRE_PtrToSolverFcn)solver->factory->methods.HYPRE_ParCSRDiagScaleSetup, 
                                                          solver->solver);
  if (SOLVER_ERROR_OCCURRED(solver))
  {
    polymec_free(solver);
    return NULL;
  }

  // Set up the virtual table.
  krylov_solver_vtable vtable = {.set_tolerances = hypre_solver_set_tolerances,
                                 .set_max_iterations = hypre_solver_set_max_iterations,
                                 .set_operator = hypre_solver_set_operator,
                                 .set_preconditioner = hypre_solver_set_pc,
                                 .solve = hypre_solver_solve,
                                 .dtor = hypre_solver_dtor};
  return krylov_solver_new("HYPRE Bi-CGSTAB", solver, vtable);
}

static krylov_solver_t* hypre_factory_special_solver(void* context,
                                                     MPI_Comm comm,
                                                     const char* solver_name,
                                                     string_string_unordered_map_t* options)
{
  hypre_solver_t* solver = polymec_malloc(sizeof(hypre_solver_t));
  solver->factory = context;

  if ((string_casecmp(solver_name, "boomerang") == 0) ||
      (string_casecmp(solver_name, "boomeramg") == 0))
  {
    solver->type = HYPRE_BOOMERANG;
    solver->factory->methods.HYPRE_BoomerAMGCreate(&solver->solver);
  }
  else
    return NULL;

  if (SOLVER_ERROR_OCCURRED(solver))
  {
    polymec_free(solver);
    return NULL;
  }

  // Set up the virtual table.
  krylov_solver_vtable vtable = {.set_tolerances = hypre_solver_set_tolerances,
                                 .set_max_iterations = hypre_solver_set_max_iterations,
                                 .set_operator = hypre_solver_set_operator,
                                 .set_preconditioner = hypre_solver_set_pc,
                                 .solve = hypre_solver_solve,
                                 .dtor = hypre_solver_dtor};
  return krylov_solver_new(solver_name, solver, vtable);
}

static void hypre_pc_free(void* pc)
{
  hypre_pc_t* p = pc;
  if ((p->dtor != NULL) && (p->pc != NULL))
    p->dtor(p->pc);
  polymec_free(p);
}

static krylov_pc_t* hypre_factory_pc(void* context,
                                     MPI_Comm comm,
                                     const char* pc_name,
                                     string_string_unordered_map_t* options)
{
  hypre_factory_t* factory = context;
  hypre_pc_t* pc = polymec_malloc(sizeof(hypre_pc_t));

  if ((string_casecmp(pc_name, "boomerang") == 0) ||
      (string_casecmp(pc_name, "boomeramg") == 0) || 
      (string_casecmp(pc_name, "amg") == 0))
  {
    factory->methods.HYPRE_BoomerAMGCreate(&pc->pc);
    pc->setup = factory->methods.HYPRE_BoomerAMGSetup;
    pc->solve = factory->methods.HYPRE_BoomerAMGSolve;
    pc->dtor = factory->methods.HYPRE_BoomerAMGDestroy;
    if (options != NULL)
    {
      int pos = 0;
      char *key, *value;
      while (string_string_unordered_map_next(options, &pos, &key, &value))
      {
        if ((strcmp(key, "num_functions") == 0) && string_is_number(value))
          factory->methods.HYPRE_BoomerAMGSetNumFunctions(pc->pc, atoi(value));
        // FIXME: etc etc
      }
    }
  }
  else if (string_casecmp(pc_name, "euclid") == 0)
  {
    HYPRE_MPI_Comm hypre_comm = comm;
    factory->methods.HYPRE_EuclidCreate(hypre_comm, &pc->pc);
    pc->setup = factory->methods.HYPRE_EuclidSetup;
    pc->solve = factory->methods.HYPRE_EuclidSolve;
    pc->dtor = factory->methods.HYPRE_EuclidDestroy;
    if (options != NULL)
    {
      int pos = 0;
      char *key, *value;
      while (string_string_unordered_map_next(options, &pos, &key, &value))
      {
        if ((strcmp(key, "level") == 0) && string_is_number(value))
          factory->methods.HYPRE_EuclidSetLevel(pc->pc, atoi(value));
        else if (strcmp(key, "bj") == 0)
          factory->methods.HYPRE_EuclidSetBJ(pc->pc, string_as_boolean(value));
        else if ((strcmp(key, "sparseA") == 0) && string_is_number(value))
          factory->methods.HYPRE_EuclidSetSparseA(pc->pc, atof(value));
        else if (strcmp(key, "rowScale") == 0)
          factory->methods.HYPRE_EuclidSetRowScale(pc->pc, string_as_boolean(value));
        else if ((strcmp(key, "ilut") == 0) && string_is_number(value))
          factory->methods.HYPRE_EuclidSetILUT(pc->pc, atof(value));
      }
    }
  }
  else if (string_casecmp(pc_name, "parasails") == 0)
  {
    HYPRE_MPI_Comm hypre_comm = comm;
    factory->methods.HYPRE_ParaSailsCreate(hypre_comm, &pc->pc);
    pc->setup = factory->methods.HYPRE_ParaSailsSetup;
    pc->solve = factory->methods.HYPRE_ParaSailsSolve;
    pc->dtor = factory->methods.HYPRE_ParaSailsDestroy;
    HYPRE_Int nlevel = 1;
    HYPRE_Real thresh = 0.1;
    if (options != NULL)
    {
      int pos = 0;
      char *key, *value;
      while (string_string_unordered_map_next(options, &pos, &key, &value))
      {
        if ((strcmp(key, "symmetry") == 0) && string_is_number(value))
          factory->methods.HYPRE_ParaSailsSetSym(pc->pc, atoi(value));
        else if ((strcmp(key, "thresh") == 0) && string_is_number(value))
          thresh = atof(value);
        else if ((strcmp(key, "nlevel") == 0) && string_is_number(value))
          nlevel = atoi(value);
        else if ((strcmp(key, "filter") == 0) && string_is_number(value))
          factory->methods.HYPRE_ParaSailsSetFilter(pc->pc, atof(value));
        else if ((strcmp(key, "loadbal") == 0) && string_is_number(value))
          factory->methods.HYPRE_ParaSailsSetLoadbal(pc->pc, atof(value));
        else if (strcmp(key, "reuse") == 0)
          factory->methods.HYPRE_ParaSailsSetReuse(pc->pc, string_as_boolean(value));
      }
    }
    factory->methods.HYPRE_ParaSailsSetParams(pc->pc, thresh, nlevel);
  }
  else
    return NULL;

  // Set up the virtual table.
  krylov_pc_vtable vtable = {.dtor = hypre_pc_free};
  return krylov_pc_new(pc_name, pc, vtable);
}

static void check_for_matrix_error(hypre_matrix_t* A, 
                                   const char* func_name)
{
  HYPRE_Int error = A->factory->methods.HYPRE_GetError();
  if (error != 0)
  {
    char msg[1024];
    A->factory->methods.HYPRE_DescribeError(error, msg);
    A->factory->methods.HYPRE_ClearAllErrors();
    polymec_error("%s: %s", func_name, msg);
  }
}

#define CHECK_FOR_MATRIX_ERROR(A) check_for_matrix_error(A, __func__);

static void hypre_matrix_assemble(void* context)
{
  hypre_matrix_t* A = context;
  A->factory->methods.HYPRE_IJMatrixAssemble(A->A);
  CHECK_FOR_MATRIX_ERROR(A)
}

static void hypre_matrix_set_values(void* context, size_t num_rows,
                                    size_t* num_columns, index_t* rows, index_t* columns,
                                    real_t* values)
{
  hypre_matrix_t* A = context;
  A->factory->methods.HYPRE_IJMatrixInitialize(A->A);
  CHECK_FOR_MATRIX_ERROR(A)

  size_t tot_num_values = 0;
  HYPRE_Int ncols[num_rows], rs[num_rows];
  for (size_t r = 0; r < num_rows; ++r)
  {
    tot_num_values += num_columns[r];
    ncols[r] = (HYPRE_Int)num_columns[r];
    rs[r] = (HYPRE_Int)rows[r];
  }
  HYPRE_Int cs[tot_num_values];
  for (size_t i = 0; i < tot_num_values; ++i)
    cs[i] = (HYPRE_Int)columns[i];
  A->factory->methods.HYPRE_IJMatrixSetValues(A->A, (HYPRE_Int)num_rows, ncols, rs, cs, values);
  CHECK_FOR_MATRIX_ERROR(A)
}

static void hypre_matrix_add_values(void* context, size_t num_rows,
                                    size_t* num_columns, index_t* rows, index_t* columns,
                                    real_t* values)
{
  hypre_matrix_t* A = context;
  A->factory->methods.HYPRE_IJMatrixInitialize(A->A);
  CHECK_FOR_MATRIX_ERROR(A)

  size_t tot_num_values = 0;
  HYPRE_Int ncols[num_rows], rs[num_rows];
  for (size_t r = 0; r < num_rows; ++r)
  {
    tot_num_values += num_columns[r];
    ncols[r] = (HYPRE_Int)num_columns[r];
    rs[r] = (HYPRE_Int)rows[r];
  }
  HYPRE_Int cs[tot_num_values];
  for (size_t i = 0; i < tot_num_values; ++i)
    cs[i] = (HYPRE_Int)columns[i];
  A->factory->methods.HYPRE_IJMatrixAddToValues(A->A, (HYPRE_Int)num_rows, ncols, rs, cs, values);
  CHECK_FOR_MATRIX_ERROR(A)
}

static void hypre_matrix_copy(void* context, void* copy)
{
  hypre_matrix_t* A = context;
  hypre_matrix_t* B = copy;
  ASSERT(B->ilow == A->ilow);
  ASSERT(B->ihigh == A->ihigh);
  ASSERT(B->block_size == A->block_size);
  B->factory->methods.HYPRE_IJMatrixInitialize(B->A);

  // Copy over the elements of A.
  size_t num_rows = (size_t)(A->ihigh - A->ilow + 1);
  size_t num_columns[num_rows], tot_num_columns = 0;
  index_t rows[num_rows];
  for (size_t r = 0; r < num_rows; ++r)
  {
    HYPRE_Int row = (HYPRE_Int)(A->ilow + r), ncols;
    A->factory->methods.HYPRE_IJMatrixGetRowCounts(A->A, 1, &row, &ncols);
    num_columns[r] = ncols;
    rows[r] = row;
    tot_num_columns += ncols;
  }
  index_t columns[tot_num_columns];
  real_t values[tot_num_columns];
  void* par_mat;
  A->factory->methods.HYPRE_IJMatrixGetObject(A->A, &par_mat);
  size_t col_offset = 0;
  for (size_t r = 0; r < num_rows; ++r)
  {
    HYPRE_Int ncols, *cols;
    HYPRE_Real* vals;
    HYPRE_Int row = (HYPRE_Int)(A->ilow + r);
    A->factory->methods.HYPRE_ParCSRMatrixGetRow(par_mat, row, &ncols, &cols, &vals);
    for (HYPRE_Int c = 0; c < ncols; ++c, ++col_offset)
    {
      columns[col_offset] = cols[c];
      values[col_offset] = (real_t)vals[c];
    }
    A->factory->methods.HYPRE_ParCSRMatrixRestoreRow(par_mat, row, &ncols, &cols, &vals);
  }
  ASSERT(col_offset == tot_num_columns);

  // Stuff the values in.
  hypre_matrix_set_values(copy, num_rows, num_columns, rows, columns, values);
  hypre_matrix_assemble(copy);
}

static void* hypre_matrix_clone(void* context)
{
  hypre_matrix_t* A = context;
  hypre_matrix_t* clone = polymec_malloc(sizeof(hypre_matrix_t));
  clone->comm = A->comm;
  clone->factory = A->factory;
  clone->ilow = A->ilow;
  clone->ihigh = A->ihigh;
  clone->block_size = A->block_size;
  if (A->block_sizes != NULL)
  {
    clone->block_sizes = polymec_malloc(sizeof(size_t) * (A->ihigh - A->ilow + 1));
    memcpy(clone->block_sizes, A->block_sizes, sizeof(size_t) * (A->ihigh - A->ilow + 1));
  }
  else
    clone->block_sizes = NULL;
  HYPRE_MPI_Comm hypre_comm = A->comm;
  clone->factory->methods.HYPRE_IJMatrixCreate(hypre_comm, 
                                               clone->ilow, clone->ihigh, 
                                               clone->ilow, clone->ihigh,
                                               &clone->A);
  clone->factory->methods.HYPRE_IJMatrixSetObjectType(clone->A, HYPRE_PARCSR);
  hypre_matrix_copy(A, clone);
  return clone;
}

static void hypre_matrix_scale(void* context, real_t scale_factor)
{
  hypre_matrix_t* A = context;

  // Get the information about the matrix's nonzero structure.
  size_t num_rows = (size_t)(A->ihigh - A->ilow + 1);
  size_t num_columns[num_rows], tot_num_columns = 0;
  index_t rows[num_rows];
  for (size_t r = 0; r < num_rows; ++r)
  {
    HYPRE_Int row = (HYPRE_Int)(A->ilow + r), ncols;
    A->factory->methods.HYPRE_IJMatrixGetRowCounts(A->A, 1, &row, &ncols);
    rows[r] = row;
    num_columns[r] = ncols;
    tot_num_columns += ncols; 
  }
  index_t columns[tot_num_columns];
  real_t values[tot_num_columns];
  void* par_mat;
  A->factory->methods.HYPRE_IJMatrixGetObject(A->A, &par_mat);
  size_t col_offset = 0;
  for (size_t r = 0; r < num_rows; ++r)
  {
    HYPRE_Int ncols, *cols;
    HYPRE_Real* vals;
    HYPRE_Int row = (HYPRE_Int)(A->ilow + r);
    A->factory->methods.HYPRE_ParCSRMatrixGetRow(par_mat, row, &ncols, &cols, &vals);
    for (HYPRE_Int c = 0; c < ncols; ++c, ++col_offset)
    {
      columns[col_offset] = cols[c];
      values[col_offset] = (real_t)(scale_factor * vals[c]);
    }
    A->factory->methods.HYPRE_ParCSRMatrixRestoreRow(par_mat, row, &ncols, &cols, &vals);
  }
  CHECK_FOR_MATRIX_ERROR(A)

  // Scale the values in the matrix.
  hypre_matrix_set_values(context, num_rows, num_columns, rows, columns, values);
  hypre_matrix_assemble(context);
}

static void hypre_vector_copy_out(void* context, real_t* local_values);
static void hypre_matrix_diag_scale(void* context, void* L, void* R)
{
  hypre_matrix_t* A = context;

  int nprocs;
  MPI_Comm_size(A->comm, &nprocs);
  if (nprocs > 1)
  {
    if (A->factory->hypre_ext == NULL)
    {
      polymec_error("hypre_matrix_diag_scale: Not supported for nprocs > 1.\n"
                    "hypre_matrix_diag_scale: Please install the HYPRE_ext library.");
    }

    // Use the HYPRE extension library to do the diagonal scaling.
    hypre_vector_t* LL = L;
    hypre_vector_t* RR = R;
    A->factory->ext_methods.HYPRE_IJMatrixDiagScale(A->A, LL->v, RR->v);
  }
  else
  {
    // We roll up our sleeves and do things manually.
    size_t num_rows = (size_t)(A->ihigh - A->ilow + 1);
    real_t Li[num_rows], Rj[num_rows];
    hypre_vector_copy_out(L, Li);
    hypre_vector_copy_out(R, Rj);

    // Get the information about the matrix's nonzero structure.
    size_t num_columns[num_rows], tot_num_columns = 0;
    index_t rows[num_rows];
    for (size_t r = 0; r < num_rows; ++r)
    {
      HYPRE_Int row = (HYPRE_Int)(A->ilow + r), ncols;
      A->factory->methods.HYPRE_IJMatrixGetRowCounts(A->A, 1, &row, &ncols);
      rows[r] = row;
      num_columns[r] = ncols;
      tot_num_columns += ncols; 
    }
    index_t columns[tot_num_columns];
    real_t values[tot_num_columns];
    void* par_mat;
    A->factory->methods.HYPRE_IJMatrixGetObject(A->A, &par_mat);
    size_t col_offset = 0;
    for (size_t r = 0; r < num_rows; ++r)
    {
      HYPRE_Int ncols, *cols;
      HYPRE_Real* vals;
      HYPRE_Int row = (HYPRE_Int)(A->ilow + r);
      A->factory->methods.HYPRE_ParCSRMatrixGetRow(par_mat, row, &ncols, &cols, &vals);
      for (HYPRE_Int c = 0; c < ncols; ++c, ++col_offset)
      {
        columns[col_offset] = cols[c];
        values[col_offset] = (real_t)(Li[r] * vals[c] * Rj[cols[c]]);
      }
      A->factory->methods.HYPRE_ParCSRMatrixRestoreRow(par_mat, row, &ncols, &cols, &vals);
    }
    CHECK_FOR_MATRIX_ERROR(A)

    // Scale the values in the matrix.
    hypre_matrix_set_values(context, num_rows, num_columns, rows, columns, values);
    hypre_matrix_assemble(context);
  }
}

static void hypre_matrix_zero(void* context)
{
  hypre_matrix_scale(context, 0.0);
}

static void hypre_matrix_add_identity(void* context, real_t scale_factor)
{
  hypre_matrix_t* A = context;

  size_t num_rows = (size_t)(A->ihigh - A->ilow + 1);
  size_t num_columns[num_rows];
  index_t rows[num_rows], columns[num_rows];
  real_t values[num_rows];
  for (size_t r = 0; r < num_rows; ++r)
  {
    num_columns[r] = 1;
    rows[r] = A->ilow + r;
    columns[r] = A->ilow + r;
    values[r] = scale_factor;
  }
  hypre_matrix_add_values(context, num_rows, num_columns, rows, columns, values); 
  hypre_matrix_assemble(context);
}

// This is used for the diagonal methods below, and defined later in the file.
static void hypre_vector_get_values(void* context, size_t num_values,
                                    index_t* indices, real_t* values);

static void hypre_matrix_set_diagonal(void* context, void* D)
{
  hypre_matrix_t* A = context;

  size_t num_rows = (size_t)(A->ihigh - A->ilow + 1);
  size_t num_columns[num_rows];
  index_t rows[num_rows], columns[num_rows];
  real_t values[num_rows];
  for (int r = 0; r < num_rows; ++r)
  {
    num_columns[r] = 1;
    rows[r] = A->ilow + r;
    columns[r] = A->ilow + r;
  }
  hypre_vector_get_values(D, num_rows, rows, values); // get values from D
  hypre_matrix_set_values(context, num_rows, num_columns, rows, columns, values); // set diag A.
  hypre_matrix_assemble(context);
}

static void hypre_matrix_add_diagonal(void* context, void* D)
{
  hypre_matrix_t* A = context;

  size_t num_rows = (size_t)(A->ihigh - A->ilow + 1);
  size_t num_columns[num_rows]; 
  index_t rows[num_rows], columns[num_rows];
  real_t values[num_rows];
  for (int r = 0; r < num_rows; ++r)
  {
    num_columns[r] = 1;
    rows[r] = A->ilow + r;
    columns[r] = A->ilow + r;
  }
  hypre_vector_get_values(D, num_rows, rows, values); // get values from D
  hypre_matrix_add_values(context, num_rows, num_columns, rows, columns, values); // add to diag A.
  hypre_matrix_assemble(context);
}

static void hypre_matrix_matvec(void* context, void* X, bool transpose, void* Y)
{
  hypre_matrix_t* A = context;
  hypre_vector_t* x = X;
  hypre_vector_t* y = Y;
  void* par_A;
  A->factory->methods.HYPRE_IJMatrixGetObject(A->A, &par_A);
  void* par_X;
  x->factory->methods.HYPRE_IJVectorGetObject(x->v, &par_X);
  void* par_Y;
  y->factory->methods.HYPRE_IJVectorGetObject(y->v, &par_Y);
  if (transpose)
    A->factory->methods.HYPRE_ParCSRMatrixMatvecT(1.0, par_A, par_X, 0.0, par_Y);
  else
    A->factory->methods.HYPRE_ParCSRMatrixMatvec(1.0, par_A, par_X, 0.0, par_Y);
}

static void hypre_matrix_get_values(void* context, size_t num_rows,
                                    size_t* num_columns, index_t* rows, index_t* columns,
                                    real_t* values)
{
  hypre_matrix_t* A = context;

  size_t tot_num_values = 0;
  HYPRE_Int ncols[num_rows], rs[num_rows];
  for (size_t r = 0; r < num_rows; ++r)
  {
    tot_num_values += num_columns[r];
    ncols[r] = (HYPRE_Int)num_columns[r];
    rs[r] = (HYPRE_Int)rows[r];
  }
  HYPRE_Int cs[tot_num_values];
  for (size_t i = 0; i < tot_num_values; ++i)
    cs[i] = (HYPRE_Int)columns[i];
  A->factory->methods.HYPRE_IJMatrixGetValues(A->A, (HYPRE_Int)num_rows, ncols, rs, cs, values);
}

static void hypre_matrix_fprintf(void* context, FILE* stream)
{
  // Write the vector to a temporary file and then write that file to 
  // the stream.
  hypre_matrix_t* A = context;
  char temp_name[FILENAME_MAX];
  FILE* temp = make_temp_file("hypre_matrixXXXXXX", temp_name);
  fclose(temp);
  A->factory->methods.HYPRE_IJMatrixPrint(A->A, temp_name);

  // Read it in as text. The actual file has a suffix with the rank.
  char actual_file[FILENAME_MAX];
  int rank;
  MPI_Comm_rank(A->comm, &rank);
  snprintf(actual_file, FILENAME_MAX-1, "%s.%05d", temp_name, rank);
  text_buffer_t* buffer = text_buffer_from_file(actual_file);
  char* text = text_buffer_to_string(buffer);
  text_buffer_free(buffer);

  // Clean up.
  remove(temp_name);
  remove(actual_file);

  // Write the text to the stream.
  fprintf(stream, "%s\n", text);
  string_free(text);
}

static void hypre_matrix_dtor(void* context)
{
  hypre_matrix_t* A = context;
  A->factory->methods.HYPRE_IJMatrixDestroy(A->A);
  A->factory = NULL;
  polymec_free(A);
}

static krylov_matrix_t* hypre_factory_matrix(void* context,
                                             matrix_sparsity_t* sparsity)
{
  hypre_matrix_t* A = polymec_malloc(sizeof(hypre_matrix_t));
  A->factory = context;
  A->comm = matrix_sparsity_comm(sparsity);
  A->block_size = 1;
  A->block_sizes = NULL;
  int rank, nprocs;
  MPI_Comm_rank(A->comm, &rank);
  MPI_Comm_size(A->comm, &nprocs);
  index_t* row_dist = matrix_sparsity_row_distribution(sparsity);
  A->ilow = (HYPRE_Int)row_dist[rank];
  A->ihigh = (HYPRE_Int)(row_dist[rank+1] - 1);
  HYPRE_MPI_Comm hypre_comm = A->comm;
  A->factory->methods.HYPRE_IJMatrixCreate(hypre_comm, 
                                           A->ilow, A->ihigh, A->ilow, A->ihigh,
                                           &A->A);
  A->factory->methods.HYPRE_IJMatrixSetObjectType(A->A, HYPRE_PARCSR);

  // Preallocate non-zero storage.
  index_t start = row_dist[rank], end = row_dist[rank+1];
  HYPRE_Int N_local = (HYPRE_Int)(end - start);
  if (nprocs == 1)
  {
    HYPRE_Int nnz[N_local];
    int rpos = 0;
    index_t row, r = 0;
    while (matrix_sparsity_next_row(sparsity, &rpos, &row))
    {
      nnz[r] = (HYPRE_Int)matrix_sparsity_num_columns(sparsity, row);
      ++r;
    }
    A->factory->methods.HYPRE_IJMatrixSetRowSizes(A->A, nnz);
  }
  else
  {
    HYPRE_Int d_nnz[N_local], o_nnz[N_local];
    memset(d_nnz, 0, sizeof(HYPRE_Int) * N_local);
    memset(o_nnz, 0, sizeof(HYPRE_Int) * N_local);
    int rpos = 0;
    index_t row, r = 0;
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
    A->factory->methods.HYPRE_IJMatrixSetDiagOffdSizes(A->A, d_nnz, o_nnz);
  }

  // We're ready to set values. 
  A->factory->methods.HYPRE_IJMatrixInitialize(A->A);

  // Set the "non-zero" values to zero initially. This constructs the specific non-zero structure.
  size_t nnz = 0;
  size_t num_columns[N_local];
  index_t rows[N_local];
  int rpos = 0;
  index_t row, r = 0;
  while (matrix_sparsity_next_row(sparsity, &rpos, &row))
  {
    rows[r] = row;
    num_columns[r] = matrix_sparsity_num_columns(sparsity, row);
    nnz += num_columns[r];
    ++r;
  }

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
  hypre_matrix_set_values(A, N_local, num_columns, rows, columns, zeros);
  hypre_matrix_assemble(A);

  // Set up the virtual table.
  krylov_matrix_vtable vtable = {.clone = hypre_matrix_clone,
                                 .copy = hypre_matrix_copy,
                                 .zero = hypre_matrix_zero,
                                 .scale = hypre_matrix_scale,
                                 .diag_scale = hypre_matrix_diag_scale,
                                 .add_identity = hypre_matrix_add_identity,
                                 .add_diagonal = hypre_matrix_add_diagonal,
                                 .set_diagonal = hypre_matrix_set_diagonal,
                                 .matvec = hypre_matrix_matvec,
                                 .set_values = hypre_matrix_set_values,
                                 .add_values = hypre_matrix_add_values,
                                 .get_values = hypre_matrix_get_values,
                                 .assemble = hypre_matrix_assemble,
                                 .fprintf = hypre_matrix_fprintf,
                                 .dtor = hypre_matrix_dtor};
  HYPRE_Int N_global = 0;
  MPI_Allreduce(&N_local, &N_global, 1, MPI_LONG_LONG, MPI_SUM, A->comm);
  return krylov_matrix_new(A, vtable, A->comm, (int)N_local, N_global);
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

static size_t matrix_block_size(void* context, index_t block_row)
{
  hypre_matrix_t* mat = context;
  if (mat->block_sizes == NULL)
    return mat->block_size;
  else return mat->block_sizes[block_row];
}

static void matrix_manipulate_fixed_blocks(void* context, size_t num_blocks,
                                           index_t* block_rows, index_t* block_columns, 
                                           real_t* block_values,
                                           void (*manipulate)(void*, size_t, size_t*, index_t*, index_t*, real_t*),
                                           bool copy_out)
{
  hypre_matrix_t* mat = context;
  size_t bs = mat->block_size;

  // We simply treat the blocks one at a time.
  for (size_t i = 0; i < num_blocks; ++i)
  {
    // Assemble the rows/columns.
    index_t block_row = block_rows[i];
    index_t block_column = block_columns[i];
    size_t num_rows = bs, num_columns[bs]; 
    index_t rows[bs], columns[bs*bs];
    {
      int l = 0;
      for (size_t j = 0; j < bs; ++j)
      {
        rows[j] = bs * block_row + j;
        num_columns[j] = bs;
        for (size_t k = 0; k < bs; ++k, ++l)
          columns[l] = bs * block_column + k;
      }
    }

    // Copy in the values if we are inserting/adding.
    real_t values[bs*bs];
    if (!copy_out)
    {
      int l = 0;
      for (size_t j = 0; j < bs; ++j)
      {
        rows[j] = bs * block_row + j;
        num_columns[j] = bs;
        for (size_t k = 0; k < bs; ++k, ++l)
          values[l] = block_values[k*bs+j];
      }
    }

    // Manipulate the values.
    manipulate(context, num_rows, num_columns, rows, columns, values);

    // Copy out the values if we are reading.
    if (copy_out)
    {
      int l = 0;
      for (size_t j = 0; j < bs; ++j)
      {
        rows[j] = bs * block_row + j;
        num_columns[j] = bs;
        for (size_t k = 0; k < bs; ++k, ++l)
          block_values[k*bs+j] = values[l];
      }
    }
  }
}

static void matrix_set_fixed_blocks(void* context, size_t num_blocks,
                                    index_t* block_rows, index_t* block_columns, 
                                    real_t* block_values)
{
  matrix_manipulate_fixed_blocks(context, num_blocks, 
                                 block_rows, block_columns, block_values,
                                 hypre_matrix_set_values, false);
}

static void matrix_add_fixed_blocks(void* context, size_t num_blocks,
                                    index_t* block_rows, index_t* block_columns, 
                                    real_t* block_values)
{
  matrix_manipulate_fixed_blocks(context, num_blocks, 
                                 block_rows, block_columns, block_values,
                                 hypre_matrix_add_values, false);
}

static void matrix_get_fixed_blocks(void* context, size_t num_blocks,
                                    index_t* block_rows, index_t* block_columns, 
                                    real_t* block_values)
{
  matrix_manipulate_fixed_blocks(context, num_blocks, 
                                 block_rows, block_columns, block_values,
                                 hypre_matrix_get_values, true);
}

static krylov_matrix_t* hypre_factory_block_matrix(void* context,
                                                   matrix_sparsity_t* sparsity,
                                                   size_t block_size)
{
  // Create a block matrix sparsity pattern with the given block size, 
  // and make a matrix from that.
  matrix_sparsity_t* block_sp = matrix_sparsity_with_block_size(sparsity, block_size);
  krylov_matrix_t* mat = hypre_factory_matrix(context, block_sp);
  matrix_sparsity_free(block_sp);

  // Set the block size and override the block matrix methods.
  hypre_matrix_t* A = mat->context;
  A->block_size = block_size;
  mat->vtable.block_size = matrix_block_size;
  mat->vtable.set_blocks = matrix_set_fixed_blocks;
  mat->vtable.add_blocks = matrix_add_fixed_blocks;
  mat->vtable.get_blocks = matrix_get_fixed_blocks;
  return mat;
}

static void matrix_manipulate_var_blocks(void* context, size_t num_blocks,
                                         index_t* block_rows, index_t* block_columns, 
                                         real_t* block_values,
                                         void (*manipulate)(void*, size_t, size_t*, index_t*, index_t*, real_t*),
                                         bool copy_out)
{
  hypre_matrix_t* mat = context;

  // We simply treat the blocks one at a time.
  for (index_t i = 0; i < num_blocks; ++i)
  {
    // Assemble the rows/columns.
    index_t block_row = block_rows[i];
    index_t block_column = block_columns[i];
    size_t bs = mat->block_sizes[block_row];
    size_t num_rows = bs;
    index_t rows[bs], columns[bs*bs];
    size_t num_columns[bs]; 
    {
      int l = 0;
      for (size_t j = 0; j < bs; ++j)
      {
        rows[j] = bs * block_row + j;
        num_columns[j] = bs;
        for (size_t k = 0; k < bs; ++k, ++l)
          columns[l] = bs * block_column + k;
      }
    }

    // Copy in the values if we are inserting/adding.
    real_t values[bs*bs];
    if (!copy_out)
    {
      int l = 0;
      for (size_t j = 0; j < bs; ++j)
      {
        rows[j] = bs * block_row + j;
        num_columns[j] = bs;
        for (size_t k = 0; k < bs; ++k, ++l)
          values[l] = block_values[k*bs+j];
      }
    }

    // Manipulate the values.
    manipulate(context, num_rows, num_columns, rows, columns, values);

    // Copy out the values if we are reading.
    if (copy_out)
    {
      int l = 0;
      for (size_t j = 0; j < bs; ++j)
      {
        rows[j] = bs * block_row + j;
        num_columns[j] = bs;
        for (size_t k = 0; k < bs; ++k, ++l)
          block_values[k*bs+j] = values[l];
      }
    }
  }
}

static void matrix_set_var_blocks(void* context, size_t num_blocks,
                                  index_t* block_rows, index_t* block_columns, 
                                  real_t* block_values)
{
  matrix_manipulate_var_blocks(context, num_blocks, 
                               block_rows, block_columns, block_values,
                               hypre_matrix_set_values, false);
}

static void matrix_add_var_blocks(void* context, size_t num_blocks,
                                  index_t* block_rows, index_t* block_columns, 
                                  real_t* block_values)
{
  matrix_manipulate_var_blocks(context, num_blocks, 
                               block_rows, block_columns, block_values,
                               hypre_matrix_add_values, false);
}

static void matrix_get_var_blocks(void* context, size_t num_blocks,
                                  index_t* block_rows, index_t* block_columns, 
                                  real_t* block_values)
{
  matrix_manipulate_var_blocks(context, num_blocks, 
                               block_rows, block_columns, block_values,
                               hypre_matrix_get_values, true);
}

static krylov_matrix_t* hypre_factory_var_block_matrix(void* context,
                                                       matrix_sparsity_t* sparsity,
                                                       size_t* block_sizes)
{
  ASSERT(block_sizes != NULL);

  // Create a block matrix sparsity pattern with the given block sizes, 
  // and make a CSR matrix from that.
  matrix_sparsity_t* block_sp = matrix_sparsity_with_block_sizes(sparsity, block_sizes);
  krylov_matrix_t* mat = hypre_factory_matrix(context, block_sp);
  matrix_sparsity_free(block_sp);

  // Set the block sizes and override the block matrix methods.
  hypre_matrix_t* A = mat->context;
  size_t N_local = (size_t)(A->ihigh - A->ilow + 1);
  A->block_sizes = polymec_malloc(sizeof(int) * N_local);
  memcpy(A->block_sizes, block_sizes, sizeof(int) * N_local);
  mat->vtable.block_size = matrix_block_size;
  mat->vtable.set_blocks = matrix_set_var_blocks;
  mat->vtable.add_blocks = matrix_add_var_blocks;
  mat->vtable.get_blocks = matrix_get_var_blocks;
  return mat;
}

static void check_for_vector_error(hypre_vector_t* v, 
                                   const char* func_name)
{
  HYPRE_Int error = v->factory->methods.HYPRE_GetError();
  if (error != 0)
  {
    char msg[1024];
    v->factory->methods.HYPRE_DescribeError(error, msg);
    v->factory->methods.HYPRE_ClearAllErrors();
    polymec_error("%s: %s", func_name, msg);
  }
}

#define CHECK_FOR_VECTOR_ERROR(v) check_for_vector_error(v, __func__);

static void hypre_vector_assemble(void* context)
{
  hypre_vector_t* v = context;
  v->factory->methods.HYPRE_IJVectorAssemble(v->v);
  CHECK_FOR_VECTOR_ERROR(v)
}

static void hypre_vector_set_values(void* context, size_t num_values,
                                    index_t* indices, real_t* values)
{
  hypre_vector_t* v = context;
  v->factory->methods.HYPRE_IJVectorInitialize(v->v);
  CHECK_FOR_VECTOR_ERROR(v)

  HYPRE_Int ind[num_values];
  for (size_t i = 0; i < num_values; ++i)
    ind[i] = (HYPRE_Int)indices[i];
  v->factory->methods.HYPRE_IJVectorSetValues(v->v, (HYPRE_Int)num_values, ind, values);
  CHECK_FOR_VECTOR_ERROR(v)
}

static void hypre_vector_add_values(void* context, size_t num_values,
                                    index_t* indices, real_t* values)
{
  hypre_vector_t* v = context;
  v->factory->methods.HYPRE_IJVectorInitialize(v->v);
  CHECK_FOR_VECTOR_ERROR(v)

  HYPRE_Int ind[num_values];
  for (size_t i = 0; i < num_values; ++i)
    ind[i] = (HYPRE_Int)indices[i];
  v->factory->methods.HYPRE_IJVectorAddToValues(v->v, (HYPRE_Int)num_values, ind, values);
  CHECK_FOR_VECTOR_ERROR(v)
}

static void hypre_vector_get_values(void* context, size_t num_values,
                                    index_t* indices, real_t* values)
{
  hypre_vector_t* v = context;

  HYPRE_Int ind[num_values];
  for (size_t i = 0; i < num_values; ++i)
    ind[i] = (HYPRE_Int)indices[i];
  v->factory->methods.HYPRE_IJVectorGetValues(v->v, (HYPRE_Int)num_values, ind, values);
}

static void hypre_vector_copy_in(void* context, real_t* local_values)
{
  hypre_vector_t* v = context;
  index_t N = v->ihigh - v->ilow + 1;
  index_t local_rows[N];
  for (index_t i = 0; i < N; ++i)
    local_rows[i] = v->ilow + i;
  hypre_vector_set_values(v, N, local_rows, local_values);
}

static void hypre_vector_copy_out(void* context, real_t* local_values)
{
  hypre_vector_t* v = context;
  index_t N = v->ihigh - v->ilow + 1;
  index_t local_rows[N];
  for (index_t i = 0; i < N; ++i)
    local_rows[i] = v->ilow + i;
  hypre_vector_get_values(v, N, local_rows, local_values);
}

static void* hypre_vector_clone(void* context)
{
  hypre_vector_t* v = context;

  // Get data from the original vector and use it to create a new one.
  hypre_vector_t* clone = polymec_malloc(sizeof(hypre_vector_t));
  clone->comm = v->comm;
  clone->factory = v->factory;
  clone->ilow = v->ilow;
  clone->ihigh = v->ihigh;
  HYPRE_MPI_Comm hypre_comm = clone->comm;
  clone->factory->methods.HYPRE_IJVectorCreate(hypre_comm, clone->ilow, clone->ihigh, &clone->v);
  clone->factory->methods.HYPRE_IJVectorSetObjectType(clone->v, HYPRE_PARCSR);
  clone->factory->methods.HYPRE_IJVectorInitialize(clone->v);

  index_t num_rows = clone->ihigh - clone->ilow + 1;
  index_t rows[num_rows];
  for (index_t r = 0; r < num_rows; ++r)
    rows[r] = clone->ilow + r;
  real_t vals[num_rows];
  hypre_vector_get_values(v, num_rows, rows, vals);
  hypre_vector_set_values(clone, num_rows, rows, vals);
  hypre_vector_assemble(clone);
  return clone;
}

static void hypre_vector_copy(void* context, void* copy)
{
  hypre_vector_t* v = context;
  hypre_vector_t* v1 = copy;

  // Get data from the original vector and use it to create a new one.
  ASSERT(v1->ilow == v->ilow);
  ASSERT(v1->ihigh = v->ihigh);
  v1->factory->methods.HYPRE_IJVectorInitialize(v1->v);

  index_t num_rows = v1->ihigh - v1->ilow + 1;
  index_t rows[num_rows];
  for (index_t r = 0; r < num_rows; ++r)
    rows[r] = v1->ilow + r;
  real_t vals[num_rows];
  hypre_vector_get_values(v, num_rows, rows, vals);
  hypre_vector_set_values(v1, num_rows, rows, vals);
  hypre_vector_assemble(v1);
}

static void hypre_vector_set_value(void* context, real_t value)
{
  hypre_vector_t* v = context;

  index_t num_rows = v->ihigh - v->ilow + 1;
  index_t rows[num_rows];
  real_t values[num_rows];
  for (index_t r = 0; r < num_rows; ++r) 
  {
    rows[r] = v->ilow + r;
    values[r] = value;
  }
  hypre_vector_set_values(context, num_rows, rows, values);
  hypre_vector_assemble(v);
}

static void hypre_vector_zero(void* context)
{
  hypre_vector_set_value(context, 0.0);
}

static void hypre_vector_scale(void* context, real_t scale_factor)
{
  hypre_vector_t* v = context;
  index_t num_rows = v->ihigh - v->ilow + 1;
  index_t rows[num_rows];
  real_t values[num_rows];
  for (index_t r = 0; r < num_rows; ++r) 
    rows[r] = v->ilow + r;
  hypre_vector_get_values(context, num_rows, rows, values);
  for (index_t r = 0; r < num_rows; ++r) 
    values[r] *= scale_factor;
  hypre_vector_set_values(context, num_rows, rows, values);
  hypre_vector_assemble(v);
}

static void hypre_vector_diag_scale(void* context, void* D)
{
  hypre_vector_t* v = context;
  index_t num_rows = v->ihigh - v->ilow + 1;
  index_t rows[num_rows];
  real_t values[num_rows], d_values[num_rows];
  for (index_t r = 0; r < num_rows; ++r) 
    rows[r] = v->ilow + r;
  hypre_vector_get_values(context, num_rows, rows, values);
  hypre_vector_get_values(D, num_rows, rows, d_values);
  for (index_t r = 0; r < num_rows; ++r) 
    values[r] *= d_values[r];
  hypre_vector_set_values(context, num_rows, rows, values);
  hypre_vector_assemble(v);
}

static real_t hypre_vector_dot(void* context, void* W)
{
  hypre_vector_t* v = context;

  // Accumulate the local part of the dot product.
  real_t local_dot = 0.0;
  index_t num_rows = v->ihigh - v->ilow + 1;
  index_t rows[num_rows];
  real_t v_values[num_rows], w_values[num_rows];
  for (index_t r = 0; r < num_rows; ++r) 
    rows[r] = v->ilow + r;
  hypre_vector_get_values(context, num_rows, rows, v_values);
  hypre_vector_get_values(W, num_rows, rows, w_values);
  for (index_t r = 0; r < num_rows; ++r) 
    local_dot += v_values[r] * w_values[r];

  // Now mash together all the parallel portions.
  real_t global_dot = 0.0;
  MPI_Allreduce(&local_dot, &global_dot, 1, MPI_REAL_T, MPI_SUM, v->comm);
  return global_dot;
}

static real_t hypre_vector_norm(void* context, int p)
{
  // SIGH. HYPRE doesn't do vector norms, so we have to do this manually.
  hypre_vector_t* v = context;
  
  // Accumulate the local part of the norm.
  real_t local_norm = 0.0;
  index_t num_rows = v->ihigh - v->ilow + 1;
  index_t rows[num_rows];
  real_t values[num_rows];
  for (index_t r = 0; r < num_rows; ++r) 
    rows[r] = v->ilow + r;
  hypre_vector_get_values(context, num_rows, rows, values);
  if (p == 0)
  {
    for (index_t r = 0; r < num_rows; ++r) 
      local_norm = MAX(local_norm, values[r]);
  }
  else if (p == 1)
  {
    for (index_t r = 0; r < num_rows; ++r) 
      local_norm += ABS(values[r]);
  }
  else
  {
    for (index_t r = 0; r < num_rows; ++r) 
      local_norm += values[r] * values[r];
  }

  // Now mash together all the parallel portions.
  real_t global_norm = 0.0;
  if (p == 0)
    MPI_Allreduce(&local_norm, &global_norm, 1, MPI_REAL_T, MPI_MAX, v->comm);
  else
    MPI_Allreduce(&local_norm, &global_norm, 1, MPI_REAL_T, MPI_SUM, v->comm);

  if (p == 2)
    global_norm = sqrt(global_norm);
  return global_norm;
}

static real_t hypre_vector_w2_norm(void* context, void* W)
{
  hypre_vector_t* v = context;
  
  // Accumulate the local part of the norm.
  real_t local_norm = 0.0;
  index_t num_rows = v->ihigh - v->ilow + 1;
  index_t rows[num_rows];
  for (index_t r = 0; r < num_rows; ++r) 
    rows[r] = v->ilow + r;
  real_t v_values[num_rows], w_values[num_rows];
  hypre_vector_get_values(context, num_rows, rows, v_values);
  hypre_vector_get_values(W, num_rows, rows, w_values);
  for (index_t i = 0; i < num_rows; ++i) 
  {
    real_t wi = w_values[i];
    real_t vi = v_values[i];
    local_norm += wi*wi*vi*vi;
  }

  // Now mash together all the parallel portions.
  real_t global_norm;
  MPI_Allreduce(&local_norm, &global_norm, 1, MPI_REAL_T, MPI_SUM, v->comm);
  return sqrt(global_norm);
}

static real_t hypre_vector_wrms_norm(void* context, void* W)
{
  hypre_vector_t* v = context;
  
  // Accumulate the local part of the norm.
  real_t local_norm = 0.0;
  index_t num_rows = v->ihigh - v->ilow + 1;
  index_t rows[num_rows];
  for (index_t r = 0; r < num_rows; ++r) 
    rows[r] = v->ilow + r;
  real_t v_values[num_rows], w_values[num_rows];
  hypre_vector_get_values(context, num_rows, rows, v_values);
  hypre_vector_get_values(W, num_rows, rows, w_values);
  for (index_t i = 0; i < num_rows; ++i) 
  {
    real_t wi = w_values[i];
    real_t vi = v_values[i];
    local_norm += wi*wi*vi*vi;
  }

  // Now mash together all the parallel portions.
  real_t local_data[2] = {local_norm, (real_t)num_rows};
  real_t global_data[2];
  MPI_Allreduce(local_data, global_data, 2, MPI_REAL_T, MPI_SUM, v->comm);
  return sqrt(global_data[0]/global_data[1]);
}

static void hypre_vector_fprintf(void* context, FILE* stream)
{
  // Write the vector to a temporary file and then write that file to 
  // the stream.
  hypre_vector_t* v = context;
  char temp_name[FILENAME_MAX];
  FILE* temp = make_temp_file("hypre_vectorXXXXXX", temp_name);
  fclose(temp);
  v->factory->methods.HYPRE_IJVectorPrint(v->v, temp_name);

  // Read it in as text. The actual file has a suffix with the rank.
  char actual_file[FILENAME_MAX];
  int rank;
  MPI_Comm_rank(v->comm, &rank);
  snprintf(actual_file, FILENAME_MAX-1, "%s.%05d", temp_name, rank);
  text_buffer_t* buffer = text_buffer_from_file(actual_file);
  char* text = text_buffer_to_string(buffer);
  text_buffer_free(buffer);

  // Clean up.
  remove(temp_name);
  remove(actual_file);

  // Write the text to the stream.
  fprintf(stream, "%s\n", text);
  string_free(text);
}

static void hypre_vector_dtor(void* context)
{
  hypre_vector_t* v = context;
  v->factory->methods.HYPRE_IJVectorDestroy(v->v);
  v->factory = NULL;
  polymec_free(v);
}

static krylov_vector_t* hypre_factory_vector(void* context,
                                             MPI_Comm comm,
                                             index_t* row_dist)
{
  hypre_vector_t* v = polymec_malloc(sizeof(hypre_vector_t));
  v->factory = context;
  v->comm = comm;
  int rank, nprocs;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &nprocs);
  v->ilow = (HYPRE_Int)row_dist[rank];
  v->ihigh = (HYPRE_Int)(row_dist[rank+1]-1);
  HYPRE_MPI_Comm hypre_comm = v->comm;
  v->factory->methods.HYPRE_IJVectorCreate(hypre_comm, v->ilow, v->ihigh, &v->v);
  v->factory->methods.HYPRE_IJVectorSetObjectType(v->v, HYPRE_PARCSR);

  // Make sure the vector is initialized.
  v->factory->methods.HYPRE_IJVectorInitialize(v->v);
  CHECK_FOR_VECTOR_ERROR(v)

  // Set all the entries of the vector to zero.
  hypre_vector_zero(v);

  // Set up the virtual table.
  krylov_vector_vtable vtable = {.clone = hypre_vector_clone,
                                 .copy = hypre_vector_copy,
                                 .zero = hypre_vector_zero,
                                 .set_value = hypre_vector_set_value,
                                 .scale = hypre_vector_scale,
                                 .diag_scale = hypre_vector_diag_scale,
                                 .set_values = hypre_vector_set_values,
                                 .add_values = hypre_vector_add_values,
                                 .get_values = hypre_vector_get_values,
                                 .copy_in = hypre_vector_copy_in,
                                 .copy_out = hypre_vector_copy_out,
                                 .assemble = hypre_vector_assemble,
                                 .dot = hypre_vector_dot,
                                 .norm = hypre_vector_norm,
                                 .w2_norm = hypre_vector_w2_norm,
                                 .wrms_norm = hypre_vector_wrms_norm,
                                 .fprintf = hypre_vector_fprintf,
                                 .dtor = hypre_vector_dtor};
  return krylov_vector_new(v, vtable, (int)(row_dist[rank+1]-row_dist[rank]), row_dist[nprocs]);
}

static void hypre_factory_dtor(void* context)
{
  hypre_factory_t* factory = context;
  if (factory->hypre_ext != NULL)
  {
    log_debug("hypre_krylov_factory: Closing HYPRE_ext library.");
    dlclose(factory->hypre_ext);
  }
  log_debug("hypre_krylov_factory: Closing HYPRE library.");
  ptr_array_free(factory->hypre_libs);
  string_array_free(factory->hypre_lib_names);
  polymec_free(factory);
}

// This is used to clean up libraries when the solver is destroyed.
static void close_lib(void* lib)
{
  dlclose(lib);
}

// Use this to retrieve symbols from dynamically loaded libraries.
#define FETCH_SYMBOL(dylib, symbol_name, function_ptr) \
  { \
    void* ptr = dlsym(dylib, symbol_name); \
    *((void**)&(function_ptr)) = ptr; \
  }

#define STR(s) SSTR(s)
#define SSTR(s) #s
krylov_factory_t* HypreFactory(const char* hypre_dir);
krylov_factory_t* HypreFactory(const char* hypre_dir)
{
  ASSERT(directory_exists(hypre_dir)); // checked by caller
  hypre_factory_t* factory = polymec_calloc(1, sizeof(hypre_factory_t));
  factory->hypre_libs = ptr_array_new();
  factory->hypre_lib_names = string_array_new();
  factory->hypre_ext = NULL;

  // Make a list of all the HYPRE libraries in our given directory.
  log_debug(STR(HypreFactory) ": Looking for HYPRE libraries in %s.", hypre_dir);
  string_array_t* all_hypre_libs = string_array_new();
  {
    string_slist_t* files_in_dir = files_within_directory(hypre_dir);
    string_slist_node_t* node = files_in_dir->front;
    while (node != NULL)
    {
      if (string_contains(node->value, "libHYPRE") && 
          string_contains(node->value, SHARED_LIBRARY_SUFFIX))
      {
        char lib_name[FILENAME_MAX+1];
        snprintf(lib_name, FILENAME_MAX, "%s/%s", hypre_dir, node->value);
        string_array_append_with_dtor(all_hypre_libs, string_dup(lib_name), string_free);
      }
      node = node->next;
    }
    string_slist_free(files_in_dir);
  }
  if (all_hypre_libs->size == 0)
  {
    log_urgent(STR(HypreFactory) ": Didn't find any HYPRE libraries!");
    goto failure;
  }

  log_debug(STR(HypreFactory) ": Found %d HYPRE libraries.", (int)all_hypre_libs->size);

  // Get the symbols.
#define FETCH_HYPRE_SYMBOL(symbol_name) \
  { \
    for (size_t ilib = 0; ilib < factory->hypre_libs->size; ++ilib) \
    { \
      FETCH_SYMBOL(factory->hypre_libs->data[ilib], #symbol_name, factory->methods.symbol_name); \
      if (factory->methods.symbol_name != NULL) \
        break; \
    } \
    if (factory->methods.symbol_name == NULL) \
    { \
      for (size_t ilib = 0; ilib < all_hypre_libs->size; ++ilib) \
      { \
        for (size_t jlib = 0; jlib < factory->hypre_libs->size; ++jlib) \
        { \
          if (strcmp(all_hypre_libs->data[ilib], factory->hypre_libs->data[jlib]) == 0) \
            continue; \
        } \
        log_debug(STR(HypreFactory) ": Looking for symbols in %s...", all_hypre_libs->data[ilib]); \
        void* lib = dlopen(all_hypre_libs->data[ilib], RTLD_NOW); \
        if (lib != NULL) \
        { \
          FETCH_SYMBOL(lib, #symbol_name, factory->methods.symbol_name); \
          if (factory->methods.symbol_name != NULL) \
          { \
            string_array_append_with_dtor(factory->hypre_lib_names, string_dup(all_hypre_libs->data[ilib]), string_free); \
            ptr_array_append_with_dtor(factory->hypre_libs, lib, close_lib); \
            break; \
          } \
          else \
            dlclose(lib); \
        } \
        else \
        { \
          char* msg = dlerror(); \
          log_debug(STR(HypreFactory) ": Couldn't open %s: %s", all_hypre_libs->data[ilib], msg); \
        } \
      } \
      if (factory->methods.symbol_name == NULL) \
      { \
        log_urgent(STR(HypreFactory) ": Couldn't find %s in any HYPRE library.", #symbol_name); \
        goto failure; \
      } \
    } \
  } \

  FETCH_HYPRE_SYMBOL(HYPRE_GetError);
  FETCH_HYPRE_SYMBOL(HYPRE_DescribeError);
  FETCH_HYPRE_SYMBOL(HYPRE_ClearAllErrors);

  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRPCGCreate);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRPCGDestroy);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRPCGSetup);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRPCGSolve);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRPCGSetTol);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRPCGSetAbsoluteTol);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRPCGSetMaxIter);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRPCGSetStopCrit);
  FETCH_HYPRE_SYMBOL(HYPRE_PCGSetTwoNorm);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRPCGSetPrecond);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRPCGGetPrecond);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRPCGGetNumIterations);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRPCGGetFinalRelativeResidualNorm);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRPCGSetPrintLevel);

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
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRGMRESSetPrintLevel);

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
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRBiCGSTABSetPrintLevel);

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
  FETCH_HYPRE_SYMBOL(HYPRE_EuclidSolve);
  FETCH_HYPRE_SYMBOL(HYPRE_EuclidSetLevel);
  FETCH_HYPRE_SYMBOL(HYPRE_EuclidSetBJ);
  FETCH_HYPRE_SYMBOL(HYPRE_EuclidSetSparseA);
  FETCH_HYPRE_SYMBOL(HYPRE_EuclidSetRowScale);
  FETCH_HYPRE_SYMBOL(HYPRE_EuclidSetILUT);

  FETCH_HYPRE_SYMBOL(HYPRE_ParaSailsCreate);
  FETCH_HYPRE_SYMBOL(HYPRE_ParaSailsDestroy);
  FETCH_HYPRE_SYMBOL(HYPRE_ParaSailsSetup);
  FETCH_HYPRE_SYMBOL(HYPRE_ParaSailsSolve);
  FETCH_HYPRE_SYMBOL(HYPRE_ParaSailsSetParams);
  FETCH_HYPRE_SYMBOL(HYPRE_ParaSailsSetFilter);
  FETCH_HYPRE_SYMBOL(HYPRE_ParaSailsSetSym);
  FETCH_HYPRE_SYMBOL(HYPRE_ParaSailsSetLoadbal);
  FETCH_HYPRE_SYMBOL(HYPRE_ParaSailsSetReuse);

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
  FETCH_HYPRE_SYMBOL(HYPRE_IJMatrixGetRowCounts);
  FETCH_HYPRE_SYMBOL(HYPRE_IJMatrixPrint);
  FETCH_HYPRE_SYMBOL(HYPRE_IJMatrixSetPrintLevel);

  FETCH_HYPRE_SYMBOL(HYPRE_IJVectorCreate);
  FETCH_HYPRE_SYMBOL(HYPRE_IJVectorDestroy);
  FETCH_HYPRE_SYMBOL(HYPRE_IJVectorInitialize);
  FETCH_HYPRE_SYMBOL(HYPRE_IJVectorSetValues);
  FETCH_HYPRE_SYMBOL(HYPRE_IJVectorAddToValues);
  FETCH_HYPRE_SYMBOL(HYPRE_IJVectorAssemble);
  FETCH_HYPRE_SYMBOL(HYPRE_IJVectorGetValues);
  FETCH_HYPRE_SYMBOL(HYPRE_IJVectorSetObjectType);
  FETCH_HYPRE_SYMBOL(HYPRE_IJVectorGetObject);
  FETCH_HYPRE_SYMBOL(HYPRE_IJVectorPrint);

  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRDiagScaleSetup);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRDiagScale);

  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRMatrixGetRow);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRMatrixRestoreRow);

  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRMatrixMatvec);
  FETCH_HYPRE_SYMBOL(HYPRE_ParCSRMatrixMatvecT);
#undef FETCH_HYPRE_SYMBOL
  log_debug(STR(HypreFactory) ": Got HYPRE symbols in %d libraries.", (int)factory->hypre_libs->size);

#if 0
  // Now try to find HYPRE_ext, the HYPRE "extension" library.
  char hypre_ext_path[FILENAME_MAX+1];
  snprintf(hypre_ext_path, FILENAME_MAX, "%s/libHYPRE_ext%s", hypre_library, SHARED_LIBRARY_SUFFIX);

  // Try to open libHYPRE_ext and mine it for symbols.
  log_debug(STR(HypreFactory) ": Opening HYPRE_ext library at %s.", hypre_ext_path);
  hypre_ext = dlopen(hypre_ext_path, RTLD_NOW);
  if (hypre_ext == NULL)
  {
    char* msg = dlerror();
    log_urgent("hypre_krylov_factory: Could not load HYPRE_ext library: ", msg);
    log_urgent("hypre_krylov_factory: Continuing without additional functionality.");
    factory->hypre_ext = NULL;
  }

#define FETCH_HYPRE_EXT_SYMBOL(symbol_name) \
  FETCH_SYMBOL(hypre_ext, #symbol_name, factory->ext_methods.symbol_name, failure);

  // Get the symbols.
  if (hypre_ext != NULL)
  {
    FETCH_HYPRE_EXT_SYMBOL(HYPRE_IJMatrixDiagScale);
    log_debug("hypre_krylov_factory: Got HYPRE_ext symbols.");
    factory->hypre_ext = hypre_ext; 
  }
#undef FETCH_HYPRE_SYMBOL
#else
    factory->hypre_ext = NULL;
#endif

  // Clean up.
  string_array_free(all_hypre_libs);

  // Construct the factory.
  krylov_factory_vtable vtable = {.pcg_solver = hypre_factory_pcg_solver,
                                  .gmres_solver = hypre_factory_gmres_solver,
                                  .bicgstab_solver = hypre_factory_bicgstab_solver,
                                  .special_solver = hypre_factory_special_solver,
                                  .preconditioner = hypre_factory_pc,
                                  .matrix = hypre_factory_matrix,
                                  .block_matrix = hypre_factory_block_matrix,
                                  .var_block_matrix = hypre_factory_var_block_matrix,
                                  .vector = hypre_factory_vector,
                                  .dtor = hypre_factory_dtor};
  return krylov_factory_new("2.10", factory, vtable);

failure:
  if (factory->hypre_ext != NULL)
    dlclose(factory->hypre_ext);
  string_array_free(all_hypre_libs);
  string_array_free(factory->hypre_lib_names);
  ptr_array_free(factory->hypre_libs);
  polymec_free(factory);
  return NULL;
}

