// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/krylov_solver.h"
#if POLYMEC_HAVE_MPI
#define USE_MPI 1
#endif
#include "lis.h"
#include "core/options.h"
#include "core/timer.h"
#include "core/string_utils.h"

// For some reason, this prototype didn't make it into lis.h.
extern LIS_INT lis_matrix_shift_diagonal(LIS_MATRIX A, LIS_SCALAR alpha);

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
} lis_solver_t;

typedef struct
{
  LIS_MATRIX A;
  int block_size, nr, nc, bnnz, nnz, *bptr, *ptr, *row, *col, *index, *bindex;
  LIS_SCALAR* values;
} lis_matrix_t;

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
  lis_matrix_t* A = op;
  solver->op = A->A;
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
  int err = lis_solve(solver->op, B, X, solver->solver);
  bool result = (err == 0);
  lis_solver_get_residualnorm(solver->solver, res_norm);
  lis_solver_get_iter(solver->solver, num_iters);
  return result;
}

static void lis_solver_dtor(void* context)
{
  lis_solver_t* solver = context;
  polymec_free(solver);
}

static krylov_solver_t* lis_factory_pcg_solver(void* context,
                                               MPI_Comm comm)
{
  lis_solver_t* solver = polymec_malloc(sizeof(lis_solver_t));
  lis_solver_create(&solver->solver);
  lis_solver_set_option("-i cg", solver->solver);
  lis_solver_set_option("-p jacobi", solver->solver);

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
  lis_solver_set_option("-p jacobi", solver->solver);

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
  lis_solver_set_option("-p jacobi", solver->solver);

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
  POLYMEC_NOT_IMPLEMENTED
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

static lis_matrix_t* matrix_new(int block_size)
{
  ASSERT(block_size > 0);
  lis_matrix_t* mat = polymec_malloc(sizeof(lis_matrix_t));
  mat->block_size = block_size;
  mat->nr = mat->nc = mat->bnnz = mat->nnz = 0;
  mat->bptr = mat->ptr = mat->row = mat->col = mat->index = mat->bindex = NULL;
  mat->values = NULL;
  return mat;
}

static void* matrix_clone(void* context)
{
  lis_matrix_t* mat = context;
  lis_matrix_t* clone = polymec_malloc(sizeof(lis_matrix_t));
  lis_matrix_duplicate(mat->A, &clone->A);
  lis_matrix_copy(mat->A, clone->A);
  return clone;
}

static void matrix_zero(void* context)
{
  lis_matrix_t* mat = context;
  memset(mat->values, 0, sizeof(LIS_SCALAR) * mat->nnz);
}

static void matrix_scale(void* context, real_t scale_factor)
{
  lis_matrix_t* mat = context;
  for (int i = 0; i < mat->nnz; ++i)
    mat->values[i] *= scale_factor;
}

static void matrix_add_identity(void* context, real_t scale_factor)
{
  lis_matrix_t* mat = context;
  lis_matrix_shift_diagonal(mat->A, (LIS_SCALAR)scale_factor);
}

static void matrix_set_diagonal(void* context, void* D)
{
  lis_matrix_t* mat = context;
  LIS_VECTOR d = D;
  LIS_INT start, end;
  lis_matrix_get_range(mat->A, &start, &end);
  for (LIS_INT i = start; i < end; ++i)
  {
    LIS_SCALAR Di;
    lis_vector_get_value(d, i, &Di);
    lis_matrix_set_value(LIS_INS_VALUE, i, i, Di, mat->A);
  }
}

static void matrix_add_diagonal(void* context, void* D)
{
  lis_matrix_t* mat = context;
  LIS_VECTOR d = D;
  LIS_INT start, end;
  lis_matrix_get_range(mat->A, &start, &end);
  for (LIS_INT i = start; i < end; ++i)
  {
    LIS_SCALAR Di;
    lis_vector_get_value(d, i, &Di);
    lis_matrix_set_value(LIS_ADD_VALUE, i, i, Di, mat->A);
  }
}

static void matrix_set_values(void* context, index_t num_rows,
                              index_t* num_columns, index_t* rows, index_t* columns,
                              real_t* values)
{
  lis_matrix_t* mat = context;
  int k = 0;
  for (int i = 0; i < num_rows; ++i)
  {
    int row = rows[i];
    for (int j = 0; j < num_columns[i]; ++j, ++k)
      lis_matrix_set_value(LIS_INS_VALUE, row, columns[k], values[k], mat->A);
  }
}

static void matrix_add_values(void* context, index_t num_rows,
                              index_t* num_columns, index_t* rows, index_t* columns,
                              real_t* values)
{
  lis_matrix_t* mat = context;
  int k = 0;
  for (int i = 0; i < num_rows; ++i)
  {
    int row = rows[i];
    for (int j = 0; j < num_columns[i]; ++j, ++k)
      lis_matrix_set_value(LIS_ADD_VALUE, row, columns[k], values[k], mat->A);
  }
}

static void matrix_start_assembly(void* context)
{
  lis_matrix_t* mat = context;
  lis_matrix_assemble(mat->A);
}

static void matrix_finish_assembly(void* context)
{
}

static void matrix_get_values(void* context, index_t num_rows,
                              index_t* num_columns, index_t* rows, index_t* columns,
                              real_t* values)
{
  POLYMEC_NOT_IMPLEMENTED
}

static void matrix_dtor(void* context)
{
  lis_matrix_t* mat = context;
  lis_matrix_destroy(mat->A);
  polymec_free(mat);
}

static krylov_matrix_t* lis_factory_matrix(void* context,
                                           adj_graph_t* sparsity)
{
  lis_factory_t* factory = context;
  MPI_Comm comm = adj_graph_comm(sparsity);
  int N_local = adj_graph_num_vertices(sparsity);
  index_t* vtx_dist = adj_graph_vertex_dist(sparsity);
  int nprocs;
  MPI_Comm_size(comm, &nprocs);
  int N_global = (int)vtx_dist[nprocs];
  int* edge_offsets = adj_graph_edge_offsets(sparsity);

  lis_matrix_t* mat = matrix_new(1);
  lis_matrix_create(comm, &mat->A);
  lis_matrix_set_size(mat->A, N_local, 0);
  mat->nnz = edge_offsets[N_local];
  mat->ptr = polymec_malloc(sizeof(LIS_INT) * (N_local+1));
  mat->index = polymec_malloc(sizeof(LIS_INT) * mat->nnz);
  mat->values = polymec_malloc(sizeof(LIS_SCALAR) * mat->nnz);
  for (int i = 0; i < N_local; ++i)
  {
    mat->ptr[i] = mat->nnz;
    mat->index[mat->nnz] = i;
    mat->nnz += 1;

    int ne = adj_graph_num_edges(sparsity, i);
    int* edges = adj_graph_edges(sparsity, i);
    for (int j = 0; j < ne; ++j, ++mat->nnz)
      mat->index[mat->nnz] = edges[j];
  }
  mat->ptr[N_local] = mat->nnz;
  memset(mat->values, 0, sizeof(LIS_SCALAR) * mat->nnz);
  lis_matrix_set_csr(mat->nnz, mat->ptr, mat->index, mat->values, mat->A);
  lis_matrix_assemble(mat->A);

  // Set up the virtual table.
  krylov_matrix_vtable vtable = {.clone = matrix_clone,
                                 .zero = matrix_zero,
                                 .scale = matrix_scale,
                                 .add_identity = matrix_add_identity,
                                 .add_diagonal = matrix_add_diagonal,
                                 .set_diagonal = matrix_set_diagonal,
                                 .set_values = matrix_set_values,
                                 .add_values = matrix_add_values,
                                 .start_assembly = matrix_start_assembly,
                                 .finish_assembly = matrix_finish_assembly,
                                 .get_values = matrix_get_values,
                                 .dtor = matrix_dtor};
  return krylov_matrix_new(mat, vtable, N_local, N_global);
}

static krylov_matrix_t* lis_factory_block_matrix(void* context,
                                                 adj_graph_t* sparsity,
                                                 int block_size)
{
  lis_factory_t* factory = context;
  MPI_Comm comm = adj_graph_comm(sparsity);
  int N_local = adj_graph_num_vertices(sparsity);
  index_t* vtx_dist = adj_graph_vertex_dist(sparsity);
  int nprocs;
  MPI_Comm_size(comm, &nprocs);
  int N_global = (int)vtx_dist[nprocs];

  // Initialize a Block Sparse Row (BSR) matrix with the given block size.
  lis_matrix_t* mat = matrix_new(block_size);
  lis_matrix_create(comm, &mat->A);
  lis_matrix_set_size(mat->A, N_local, 0);
  int* adj = adj_graph_adjacency(sparsity);
  mat->bnnz = adj[N_local];
  mat->nnz = block_size * block_size * mat->bnnz;
  mat->nr = N_local, mat->nc = N_global;
  mat->bptr = malloc(sizeof(LIS_INT) * (mat->nr+1));
  mat->bindex = malloc(sizeof(LIS_INT) * mat->bnnz);
  mat->values = malloc(sizeof(LIS_SCALAR) * mat->nnz);
  int boffset = 0;
  for (int i = 0; i < mat->nr; ++i)
  {
    mat->bptr[i] = i;
    mat->bindex[boffset] = i;
    ++boffset;

    int ne = adj_graph_num_edges(sparsity, i);
    int* edges = adj_graph_edges(sparsity, i);
    for (int j = 0; j < ne; ++j, ++boffset)
    {
      int jj = edges[j];
      mat->bindex[boffset] = jj;
      memset(&mat->values[block_size*block_size*boffset], 0, sizeof(LIS_SCALAR) * block_size * block_size);
    }
  }
  mat->bptr[mat->nr] = mat->nr;
  lis_matrix_set_bsr(block_size, block_size, mat->bnnz, mat->bptr, mat->bindex, mat->values, mat->A);
  lis_matrix_assemble(mat->A);

  // Set up the virtual table.
  krylov_matrix_vtable vtable = {.clone = matrix_clone,
                                 .zero = matrix_zero,
                                 .scale = matrix_scale,
                                 .add_identity = matrix_add_identity,
                                 .add_diagonal = matrix_add_diagonal,
                                 .set_diagonal = matrix_set_diagonal,
                                 .set_values = matrix_set_values,
                                 .add_values = matrix_add_values,
                                 .start_assembly = matrix_start_assembly,
                                 .finish_assembly = matrix_finish_assembly,
                                 .get_values = matrix_get_values,
                                 .dtor = matrix_dtor};
  return krylov_matrix_new(mat, vtable, N_local, N_global);
}

static krylov_matrix_t* lis_factory_var_block_matrix(void* context,
                                                     adj_graph_t* sparsity,
                                                     int* block_sizes)
{
  lis_factory_t* factory = context;
  MPI_Comm comm = adj_graph_comm(sparsity); 
  int N_local = adj_graph_num_vertices(sparsity);
  index_t* vtx_dist = adj_graph_vertex_dist(sparsity);
  int nprocs;
  MPI_Comm_size(comm, &nprocs);
  int N_global = (int)vtx_dist[nprocs];

  lis_matrix_t* mat = matrix_new(1);
  lis_matrix_create(comm, &mat->A);
  lis_matrix_set_size(mat->A, N_local, 0);
  mat->nr = N_local; 
  mat->nc = N_global;
  int* adj = adj_graph_adjacency(sparsity);
  mat->bnnz = adj[N_local];
  mat->row = malloc(sizeof(LIS_INT) * (mat->nr+1));
  mat->col = malloc(sizeof(LIS_INT) * (mat->nc+1));
  mat->ptr = malloc(sizeof(LIS_INT) * (mat->bnnz+1));
  mat->bindex = malloc(sizeof(LIS_INT) * mat->bnnz);
  mat->bptr = malloc(sizeof(LIS_INT) * (mat->nr+1));
  LIS_INT boffset = 0, row_offset = 0;
  mat->nnz = 0;
  for (int i = 0; i < N_local; ++i)
  {
    int block_size = block_sizes[i];
    mat->row[i] = row_offset;
    mat->col[i] = row_offset;
    row_offset += block_size;
    mat->bptr[i] = i;
    mat->nnz += block_size * block_size;

    int ne = adj_graph_num_edges(sparsity, i);
    int* edges = adj_graph_edges(sparsity, i);
    for (int j = 0; j < ne; ++j, ++boffset)
    {
      int jj = edges[j];
      mat->bindex[boffset] = jj;
      mat->ptr[boffset] = mat->nnz;
      mat->nnz += block_size * block_size;
    }
  }
  mat->ptr[mat->bnnz] = mat->nnz;
  mat->values = malloc(sizeof(LIS_SCALAR) * mat->nnz);
  memset(mat->values, 0, sizeof(LIS_SCALAR) * mat->nnz);
  lis_matrix_set_vbr(mat->nnz, mat->nr, mat->nc, mat->bnnz, mat->row, mat->col, mat->ptr, mat->bptr, mat->bindex, mat->values, mat->A);
  lis_matrix_assemble(mat->A);

  // Set up the virtual table.
  krylov_matrix_vtable vtable = {.clone = matrix_clone,
                                 .zero = matrix_zero,
                                 .scale = matrix_scale,
                                 .add_identity = matrix_add_identity,
                                 .add_diagonal = matrix_add_diagonal,
                                 .set_diagonal = matrix_set_diagonal,
                                 .set_values = matrix_set_values,
                                 .add_values = matrix_add_values,
                                 .start_assembly = matrix_start_assembly,
                                 .finish_assembly = matrix_finish_assembly,
                                 .get_values = matrix_get_values,
                                 .dtor = matrix_dtor};
  return krylov_matrix_new(mat, vtable, N_local, N_global);
}

static void* vector_clone(void* context)
{
  LIS_VECTOR v = context;
  LIS_VECTOR clone;
  lis_vector_duplicate(v, &clone);
  lis_vector_copy(v, clone);
  return clone;
}

static void vector_zero(void* context)
{
  LIS_VECTOR v = context;
  lis_vector_set_all(0.0, v);
}

static void vector_set_value(void* context, real_t value)
{
  LIS_VECTOR v = context;
  lis_vector_set_all((LIS_SCALAR)value, v);
}

static void vector_scale(void* context, real_t scale_factor)
{
  LIS_VECTOR v = context;
  lis_vector_scale((LIS_SCALAR)scale_factor, v);
}

static void vector_set_values(void* context, index_t num_values,
                              index_t* indices, real_t* values)
{
  LIS_VECTOR v = context;
  lis_vector_set_values(LIS_INS_VALUE, (LIS_INT)num_values, 
                        (LIS_INT*)indices, values, v);
}

static void vector_add_values(void* context, index_t num_values,
                              index_t* indices, real_t* values)
{
  LIS_VECTOR v = context;
  lis_vector_set_values(LIS_ADD_VALUE, (LIS_INT)num_values, 
                        (LIS_INT*)indices, values, v);
}

static void vector_start_assembly(void* context)
{
}

static void vector_finish_assembly(void* context)
{
}

static void vector_get_values(void* context, index_t num_values,
                              index_t* indices, real_t* values)
{
  LIS_VECTOR v = context;
  for (index_t i = 0; i < num_values; ++i)
    lis_vector_get_values(v, indices[i], 1, &values[i]);
}

static real_t vector_norm(void* context, int p)
{
  LIS_VECTOR v = context;
  real_t norm;
  if (p == 0)
    lis_vector_nrmi(v, &norm);
  else if (p == 1)
    lis_vector_nrm1(v, &norm);
  else // (p == 2)
    lis_vector_nrm2(v, &norm);
  return norm;
}

static void vector_dtor(void* context)
{
  LIS_VECTOR v = context;
  lis_vector_destroy(v);
}

static krylov_vector_t* lis_factory_vector(void* context,
                                           adj_graph_t* dist_graph)
{
  lis_factory_t* factory = context;
  int N_local = adj_graph_num_vertices(dist_graph);
  index_t* vtx_dist = adj_graph_vertex_dist(dist_graph);
  MPI_Comm comm = adj_graph_comm(dist_graph);
  int nprocs, rank;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_size(comm, &rank);
  int N_global = (int)vtx_dist[nprocs];

  LIS_VECTOR v;
  lis_vector_create(comm, &v);
  lis_vector_set_size(v, N_local, N_global);

  // Set up the virtual table.
  krylov_vector_vtable vtable = {.clone = vector_clone,
                                 .zero = vector_zero,
                                 .set_value = vector_set_value,
                                 .scale = vector_scale,
                                 .set_values = vector_set_values,
                                 .add_values = vector_add_values,
                                 .start_assembly = vector_start_assembly,
                                 .finish_assembly = vector_finish_assembly,
                                 .get_values = vector_get_values,
                                 .norm = vector_norm,
                                 .dtor = vector_dtor};
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
                                  .var_block_matrix = lis_factory_var_block_matrix,
                                  .vector = lis_factory_vector,
                                  .dtor = lis_factory_dtor};
  return krylov_factory_new("LIS", factory, vtable);
}

