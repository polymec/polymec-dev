// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "lis.h"
#include "core/krylov_solver.h"
#include "core/options.h"
#include "core/timer.h"
#include "core/string_utils.h"

//------------------------------------------------------------------------
// This file implements the LIS Krylov solver.
//------------------------------------------------------------------------

typedef struct
{
  MPI_Comm comm;
} lis_factory_t;

typedef struct
{
  LIS_SOLVER solver;
  LIS_MATRIX op;
  LIS_PRECON pc;
} lis_solver_t;

static void lis_solver_set_tolerances(void* context,
                                      real_t rel_tol,
                                      real_t abs_tol,
                                      real_t div_tol)

{
  lis_solver_t* solver = context;
  char command[129];
  snprintf(command, 128, "-conv_cond 0 -tol %g", rel_tol);
  lis_solver_set_option(command, solver->solver);
}

static void lis_solver_set_max_iterations(void* context,
                                          int max_iters)
{
  lis_solver_t* solver = context;
  char command[129];
  snprintf(command, 128, "-maxiter %d", max_iters);
  lis_solver_set_option(command, solver->solver);
}

static void lis_solver_set_operator(void* context,
                                    void* op)
{
  lis_solver_t* solver = context;
  solver->op = op;
}

static void lis_solver_set_pc(void* context,
                              void* pc)
{
  lis_solver_t* solver = context;
  char* text = pc;
  lis_solver_set_option(text, solver->solver);
  polymec_free(text);
}

static bool lis_solver_solve(void* context,
                             void* b,
                             void* x,
                             real_t* res_norm,
                             int* num_iters)
{
  lis_solver_t* solver = context;
  LIS_VECTOR B = b;
  LIS_VECTOR X = x;
  int err = lis_solve_kernel(solver->op, B, X, solver->solver, solver->pc);
  bool result = (err == 0);
  lis_solver_get_residualnorm(solver->solver, res_norm);
  lis_solver_get_iter(solver->solver, num_iters);
  return result;
}

static void lis_solver_dtor(void* context)
{
  lis_solver_t* solver = context;
  solver->factory->methods.KSPDestroy(&solver->ksp);
  solver->factory = NULL;
  polymec_free(solver);
}

static krylov_solver_t* lis_factory_pcg_solver(void* context,
                                               MPI_Comm comm)
{
  lis_solver_t* solver = polymec_malloc(sizeof(lis_solver_t));
  lis_solver_create(&solver->solver);
  lis_solver_set_option("-i cg", solver->solver);
  lis_precon_create(solver->solver, &solver->pc);
  lis_solver_set_option("-p jacobi");

  // Set up the virtual table.
  krylov_solver_vtable vtable = {.set_tolerances = lis_solver_set_tolerances,
                                 .set_max_iterations = lis_solver_set_max_iterations,
                                 .set_operator = lis_solver_set_operator,
                                 .set_preconditioner = lis_solver_set_pc,
                                 .solve = lis_solver_solve,
                                 .dtor = lis_solver_dtor};

  return krylov_solver_new("LIS PCG", solver, vtable);
}

static krylov_solver_t* lis_factory_gmres_solver(void* context,
                                                 MPI_Comm comm,
                                                 int krylov_dimension)
{
  lis_solver_t* solver = polymec_malloc(sizeof(lis_solver_t));
  lis_solver_create(&solver->solver);
  lis_solver_set_option("-i gmres", solver->solver);
  lis_precon_create(solver->solver, &solver->pc);
  lis_solver_set_option("-p jacobi");

  // Set up the virtual table.
  krylov_solver_vtable vtable = {.set_tolerances = lis_solver_set_tolerances,
                                 .set_max_iterations = lis_solver_set_max_iterations,
                                 .set_operator = lis_solver_set_operator,
                                 .set_preconditioner = lis_solver_set_pc,
                                 .solve = lis_solver_solve,
                                 .dtor = lis_solver_dtor};
  return krylov_solver_new("LIS GMRES", solver, vtable);
}

static krylov_solver_t* lis_factory_bicgstab_solver(void* context,
                                                    MPI_Comm comm)
{
  lis_solver_t* solver = polymec_malloc(sizeof(lis_solver_t));
  lis_solver_create(&solver->solver);
  lis_solver_set_option("-i bicgstab", solver->solver);
  lis_precon_create(solver->solver, &solver->pc);
  lis_solver_set_option("-p jacobi");

  // Set up the virtual table.
  krylov_solver_vtable vtable = {.set_tolerances = lis_solver_set_tolerances,
                                 .set_max_iterations = lis_solver_set_max_iterations,
                                 .set_operator = lis_solver_set_operator,
                                 .set_preconditioner = lis_solver_set_pc,
                                 .solve = lis_solver_solve,
                                 .dtor = lis_solver_dtor};
  return krylov_solver_new("LIS BiCGSTAB", solver, vtable);
}

static krylov_solver_t* lis_factory_special_solver(void* context,
                                                   MPI_Comm comm,
                                                   const char* solver_name,
                                                   string_string_unordered_map_t* options)
{
  // FIXME
  return NULL;
}

static krylov_pc_t* lis_factory_pc(void* context,
                                   MPI_Comm comm,
                                   const char* pc_name,
                                   string_string_unordered_map_t* options)
{
  // We store the preconditioner as a text string.
  char* text = polymec_malloc(sizeof(char) * 129);
  snprintf(text, 128, "-p %s", pc_name);

  // Set up the virtual table.
  krylov_pc_vtable vtable = {.dtor = NULL}; 
  return krylov_pc_new(pc_name, text, vtable);
}

static void* lis_matrix_clone(void* context)
{
  LIS_MATRIX A = context;
  LIS_MATRIX clone;
  lis_matrix_duplicate(A, &clone);
  lis_matrix_copy(A, clone);
  return clone;
}

static void lis_matrix_zero(void* context)
{
  LIS_MATRIX A = context;
  // FIXME
}

static void lis_matrix_scale(void* context, real_t scale_factor)
{
  LIS_MATRIX A = context;
  // FIXME
}

static void lis_matrix_add_identity(void* context, real_t scale_factor)
{
  LIS_MATRIX A = context;
  // FIXME
}

static void lis_matrix_set_diagonal(void* context, void* D)
{
  LIS_MATRIX A = context;
  // FIXME
}

static void lis_matrix_add_diagonal(void* context, void* D)
{
  LIS_MATRIX A = context;
  // FIXME
}

static void lis_matrix_insert_values(void* context, index_t num_rows,
                                     index_t* num_columns, index_t* rows, index_t* columns,
                                     real_t* values, InsertMode insert_mode)
{
  // FIXME
}

static void lis_matrix_set_values(void* context, index_t num_rows,
                                  index_t* num_columns, index_t* rows, index_t* columns,
                                  real_t* values)
{
  LIS_MATRIX A = context;
  // FIXME
}

static void lis_matrix_add_values(void* context, index_t num_rows,
                                  index_t* num_columns, index_t* rows, index_t* columns,
                                  real_t* values)
{
  LIS_MATRIX A = context;
  // FIXME
}

static void lis_matrix_start_assembly(void* context)
{
  LIS_MATRIX A = context;
  lis_matrix_assemble(A);
}

static void lis_matrix_finish_assembly(void* context)
{
}

static void lis_matrix_get_values(void* context, index_t num_rows,
                                  index_t* num_columns, index_t* rows, index_t* columns,
                                  real_t* values)
{
  LIS_MATRIX A = context;
  // FIXME
}

static void lis_matrix_dtor(void* context)
{
  LIS_MATRIX A = context;
  lis_matrix_destroy(A);
}

static krylov_matrix_t* lis_factory_matrix(void* context,
                                           adj_graph_t* sparsity)
{
  lis_factory_t* factory = context;
  MPI_Comm = adj_graph_comm(sparsity); 

  LIS_MATRIX A;
  lis_matrix_create(comm, &A);
  // FIXME
  int N_local = 0, N_global = 0;

  // Set up the virtual table.
  krylov_matrix_vtable vtable = {.clone = lis_matrix_clone,
                                 .zero = lis_matrix_zero,
                                 .scale = lis_matrix_scale,
                                 .add_identity = lis_matrix_add_identity,
                                 .add_diagonal = lis_matrix_add_diagonal,
                                 .set_diagonal = lis_matrix_set_diagonal,
                                 .set_values = lis_matrix_set_values,
                                 .add_values = lis_matrix_add_values,
                                 .start_assembly = lis_matrix_start_assembly,
                                 .finish_assembly = lis_matrix_finish_assembly,
                                 .get_values = lis_matrix_get_values,
                                 .dtor = lis_matrix_dtor};
  return krylov_matrix_new(A, vtable, N_local, N_global);
}

static krylov_matrix_t* lis_factory_block_matrix(void* context,
                                                 adj_graph_t* sparsity,
                                                 int block_size)
{
  lis_factory_t* factory = context;
  MPI_Comm = adj_graph_comm(sparsity); 

  LIS_MATRIX A;
  lis_matrix_create(comm, &A);
  // FIXME
  lis_matrix_set_blocksize(A, block_size, block_size, NULL, NULL);
  int N_local = 0, N_global = 0;

  // Set up the virtual table.
  krylov_matrix_vtable vtable = {.clone = lis_matrix_clone,
                                 .zero = lis_matrix_zero,
                                 .scale = lis_matrix_scale,
                                 .add_identity = lis_matrix_add_identity,
                                 .add_diagonal = lis_matrix_add_diagonal,
                                 .set_diagonal = lis_matrix_set_diagonal,
                                 .set_values = lis_matrix_set_values,
                                 .add_values = lis_matrix_add_values,
                                 .start_assembly = lis_matrix_start_assembly,
                                 .finish_assembly = lis_matrix_finish_assembly,
                                 .get_values = lis_matrix_get_values,
                                 .dtor = lis_matrix_dtor};
  return krylov_matrix_new(A, vtable, N_local, N_global);
}

static void* lis_vector_clone(void* context)
{
  LIS_VECTOR v = context;
  LIS_VECTOR clone;
  lis_matrix_duplicate(v, &clone);
  lis_matrix_copy(v, clone);
  return clone;
}

static void lis_vector_zero(void* context)
{
  LIS_VECTOR v = context;
  // FIXME
}

static void lis_vector_set_value(void* context, real_t value)
{
  LIS_VECTOR v = context;
  // FIXME
}

static void lis_vector_scale(void* context, real_t scale_factor)
{
  LIS_VECTOR v = context;
  // FIXME
}

static void lis_vector_set_values(void* context, index_t num_values,
                                  index_t* indices, real_t* values)
{
  LIS_VECTOR v = context;
  // FIXME
}

static void lis_vector_add_values(void* context, index_t num_values,
                                  index_t* indices, real_t* values)
{
  LIS_VECTOR v = context;
  // FIXME
}

static void lis_vector_start_assembly(void* context)
{
  LIS_VECTOR v = context;
  lis_vector_assemble(v);
}

static void lis_vector_finish_assembly(void* context)
{
}

static void lis_vector_get_values(void* context, index_t num_values,
                                  index_t* indices, real_t* values)
{
  LIS_VECTOR v = context;
  // FIXME
}

static real_t lis_vector_norm(void* context, int p)
{
  LIS_VECTOR v = context;
  // FIXME
  return -FLT_MAX;
}

static void lis_vector_dtor(void* context)
{
  LIS_VECTOR v = context;
  lis_vector_destroy(v);
}

static krylov_vector_t* lis_factory_vector(void* context,
                                           adj_graph_t* dist_graph)
{
  lis_factory_t* factory = context;
  MPI_Comm = adj_graph_comm(sparsity); 

  LIS_VECTOR v;
  lis_vector_create(comm, &v);
  // FIXME

  int rank, nprocs;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &nprocs);
  index_t* vtx_dist = adj_graph_vertex_dist(dist_graph);

  // Set up the virtual table.
  krylov_vector_vtable vtable = {.clone = lis_vector_clone,
                                 .zero = lis_vector_zero,
                                 .set_value = lis_vector_set_value,
                                 .scale = lis_vector_scale,
                                 .set_values = lis_vector_set_values,
                                 .add_values = lis_vector_add_values,
                                 .start_assembly = lis_vector_start_assembly,
                                 .finish_assembly = lis_vector_finish_assembly,
                                 .get_values = lis_vector_get_values,
                                 .norm = lis_vector_norm,
                                 .dtor = lis_vector_dtor};
  return krylov_vector_new(v, vtable, vtx_dist[rank+1]-vtx_dist[rank], vtx_dist[nprocs]);
}

static void lis_factory_dtor(void* context)
{
  lis_factory_t* factory = context;
  polymec_free(factory);
}

krylov_factory_t* lis_krylov_factory()
{
  lis_factory_t* factory = polymec_malloc(sizeof(lis_factory_t));

  // Construct the factory.
  krylov_factory_vtable vtable = {.pcg_solver = lis_factory_pcg_solver,
                                  .gmres_solver = lis_factory_gmres_solver,
                                  .bicgstab_solver = lis_factory_bicgstab_solver,
                                  .special_solver = lis_factory_special_solver,
                                  .preconditioner = lis_factory_pc,
                                  .matrix = lis_factory_matrix,
                                  .block_matrix = lis_factory_block_matrix,
                                  .vector = lis_factory_vector,
                                  .dtor = lis_factory_dtor};
  return krylov_factory_new("LIS", factory, vtable);
}

