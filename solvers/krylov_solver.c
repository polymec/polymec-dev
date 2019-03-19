// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <float.h>
#include "core/polymec.h"
#include "core/timer.h"
#include "core/array_utils.h"
#include "solvers/krylov_solver.h"

extern void distribute_matrix_sparsity(matrix_sparsity_t** sparsity,
                                       MPI_Comm comm,
                                       index_t* row_distribution);

//------------------------------------------------------------------------
//                Krylov data types and virtual tables
//------------------------------------------------------------------------

struct krylov_solver_t
{
  char* name;
  void* context;
  krylov_solver_vtable vtable;
  krylov_pc_t* pc;
  krylov_matrix_t* op;

  // Scaled solver stuff.
  krylov_matrix_t* scaled_op;
  krylov_vector_t* scaled_b;
  krylov_vector_t* s2_inv;
};

struct krylov_pc_t
{
  char* name;
  void* context;
  krylov_pc_vtable vtable;
};

struct krylov_matrix_t
{
  void* context;
  krylov_matrix_vtable vtable;
  MPI_Comm comm;
  size_t num_local_rows;
  size_t num_global_rows;
};

struct krylov_vector_t
{
  void* context;
  krylov_vector_vtable vtable;
  size_t local_size;
  size_t global_size;
};

struct krylov_factory_t
{
  char* name;
  void* context;
  krylov_factory_vtable vtable;
};

//------------------------------------------------------------------------
//                          Krylov solver
//------------------------------------------------------------------------

krylov_solver_t* krylov_solver_new(const char* name,
                                   void* context,
                                   krylov_solver_vtable vtable)
{
  ASSERT(vtable.set_tolerances != NULL);
  ASSERT(vtable.set_max_iterations != NULL);
  ASSERT(vtable.set_operator != NULL);
  ASSERT(vtable.set_preconditioner != NULL);
  ASSERT(vtable.solve != NULL);
  krylov_solver_t* solver = polymec_malloc(sizeof(krylov_solver_t));
  solver->name = string_dup(name);
  solver->context = context;
  solver->vtable = vtable;
  solver->pc = NULL;
  solver->op = NULL;
  solver->scaled_op = NULL;
  solver->scaled_b = NULL;
  solver->s2_inv = NULL;
  return solver;
}

void krylov_solver_free(krylov_solver_t* solver)
{
  if ((solver->context != NULL) && (solver->vtable.dtor != NULL))
    solver->vtable.dtor(solver->context);
  if (solver->pc != NULL)
    krylov_pc_free(solver->pc);
  if (solver->scaled_op != NULL)
    krylov_matrix_free(solver->scaled_op);
  if (solver->scaled_b != NULL)
    krylov_vector_free(solver->scaled_b);
  if (solver->s2_inv != NULL)
    krylov_vector_free(solver->s2_inv);
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
  START_FUNCTION_TIMER();
  solver->op = op;
  solver->vtable.set_operator(solver->context, op->context);
  STOP_FUNCTION_TIMER();
}

void krylov_solver_set_preconditioner(krylov_solver_t* solver,
                                      krylov_pc_t* preconditioner)
{
  if (solver->pc != NULL)
    krylov_pc_free(solver->pc);
  solver->pc = preconditioner;
}

krylov_pc_t* krylov_solver_preconditioner(krylov_solver_t* solver)
{
  return solver->pc;
}

bool krylov_solver_solve(krylov_solver_t* solver,
                         krylov_vector_t* b,
                         krylov_vector_t* x,
                         real_t* residual_norm,
                         int* num_iterations)
{
  START_FUNCTION_TIMER();
  ASSERT(solver->op != NULL);
  bool solved = solver->vtable.solve(solver->context, b->context, x->context,
                                     residual_norm, num_iterations);
  STOP_FUNCTION_TIMER();
  return solved;
}

bool krylov_solver_solve_scaled(krylov_solver_t* solver,
                                krylov_vector_t* b,
                                krylov_vector_t* s1,
                                krylov_vector_t* s2,
                                krylov_vector_t* x,
                                real_t* residual_norm,
                                int* num_iterations)
{
  START_FUNCTION_TIMER();
  ASSERT(solver->op != NULL);

  // If we're not using the scaling matrices, do an unscaled solve.
  bool solved;
  if ((s1 == NULL) && (s2 == NULL))
    solved = krylov_solver_solve(solver, b, x, residual_norm, num_iterations);
  else
  {
    // Make sure s2 inverse is computed if s2 is given.
    if (s2 != NULL)
    {
      if (solver->s2_inv == NULL)
        solver->s2_inv = krylov_vector_clone(s2);

      size_t N_local = krylov_vector_local_size(s2);
      real_t s2_data[N_local], s2_inv_data[N_local];
      krylov_vector_copy_out(s2, s2_data);
      for (int i = 0; i < (int)N_local; ++i)
        s2_inv_data[i] = 1.0/s2_data[i];
      krylov_vector_copy_in(solver->s2_inv, s2_inv_data);
    }
    else
    {
      if (solver->s2_inv != NULL)
      {
        krylov_vector_free(solver->s2_inv);
        solver->s2_inv = NULL;
      }
    }

    // Calculate the scaled operator matrix s1 * A * s2_inv.
    if (solver->scaled_op == NULL)
      solver->scaled_op = krylov_matrix_clone(solver->op);
    else
      krylov_matrix_copy(solver->op, solver->scaled_op);
    krylov_matrix_diag_scale(solver->scaled_op, s1, solver->s2_inv);

    // Calculate the scaled right-hand side s1 * b.
    if (solver->scaled_b == NULL)
      solver->scaled_b = krylov_vector_clone(b);
    else
      krylov_vector_copy(b, solver->scaled_b);
    if (s1 != NULL)
      krylov_vector_diag_scale(solver->scaled_b, s1);

    // Solve the scaled system (s1 * A * s2_inv) * (s2 * x) = (s1 * b).
    // The residual norm is ||s1 * P^{-1} * (b - A * x)||_2, where P is
    // the preconditioner matrix for the solver.
    solver->vtable.set_operator(solver->context, solver->scaled_op->context);
    solved = solver->vtable.solve(solver->context, solver->scaled_b->context,
                                  x->context, residual_norm, num_iterations);
    solver->vtable.set_operator(solver->context, solver->op->context);

    // Transform (s2 * x) -> x and return.
    if (solved && (solver->s2_inv != NULL))
      krylov_vector_diag_scale(x, solver->s2_inv);
  }
  STOP_FUNCTION_TIMER();
  return solved;
}

//------------------------------------------------------------------------
//                          Krylov preconditioner
//------------------------------------------------------------------------

krylov_pc_t* krylov_pc_new(const char* name,
                           void* context,
                           krylov_pc_vtable vtable)
{
  krylov_pc_t* pc = polymec_malloc(sizeof(krylov_pc_t));
  pc->name = string_dup(name);
  pc->context = context;
  pc->vtable = vtable;
  return pc;
}

void krylov_pc_free(krylov_pc_t* preconditioner)
{
  if ((preconditioner->context != NULL) && (preconditioner->vtable.dtor != NULL))
    preconditioner->vtable.dtor(preconditioner->context);
  polymec_free(preconditioner);
}

char* krylov_pc_name(krylov_solver_t* preconditioner)
{
  return preconditioner->name;
}

//------------------------------------------------------------------------
//                          Krylov matrix
//------------------------------------------------------------------------

krylov_matrix_t* krylov_matrix_new(void* context,
                                   krylov_matrix_vtable vtable,
                                   MPI_Comm comm,
                                   size_t num_local_rows,
                                   size_t num_global_rows)
{
  ASSERT(vtable.clone != NULL);
  ASSERT(vtable.copy != NULL);
  ASSERT(vtable.zero != NULL);
  ASSERT(vtable.scale != NULL);
  ASSERT(vtable.diag_scale != NULL);
  ASSERT(vtable.add_diagonal != NULL);
  ASSERT(vtable.set_diagonal != NULL);
  ASSERT(vtable.matvec != NULL);
  ASSERT(vtable.set_values != NULL);
  ASSERT(vtable.add_values != NULL);
  ASSERT(vtable.get_values != NULL);
  ASSERT(num_local_rows > 0);
  ASSERT(num_global_rows >= num_local_rows);
  krylov_matrix_t* A = polymec_malloc(sizeof(krylov_matrix_t));
  A->context = context;
  A->vtable = vtable;
  A->comm = comm;
  A->num_local_rows = num_local_rows;
  A->num_global_rows = num_global_rows;
  return A;
}

//------------------------------------------------------------------------
//                 Matrix Market file format logic
//------------------------------------------------------------------------
// This code was adapted from http://math.nist.gov/MatrixMarket.
//------------------------------------------------------------------------

#define MM_MAX_LINE_LENGTH 1025
#define MatrixMarketBanner "%%MatrixMarket"
#define MM_MAX_TOKEN_LENGTH 64

#define mm_is_matrix(typecode)  ((typecode)[0]=='M')

#define mm_is_sparse(typecode)  ((typecode)[1]=='C')
#define mm_is_coordinate(typecode)((typecode)[1]=='C')
#define mm_is_dense(typecode)  ((typecode)[1]=='A')
#define mm_is_array(typecode)  ((typecode)[1]=='A')

#define mm_is_complex(typecode)  ((typecode)[2]=='C')
#define mm_is_real(typecode)    ((typecode)[2]=='R')
#define mm_is_pattern(typecode)  ((typecode)[2]=='P')
#define mm_is_integer(typecode) ((typecode)[2]=='I')

#define mm_is_symmetric(typecode)((typecode)[3]=='S')
#define mm_is_general(typecode)  ((typecode)[3]=='G')
#define mm_is_skew(typecode)  ((typecode)[3]=='K')
#define mm_is_hermitian(typecode)((typecode)[3]=='H')

#define MM_MTX_STR    "matrix"
#define MM_ARRAY_STR  "array"
#define MM_DENSE_STR  "array"
#define MM_COORDINATE_STR "coordinate"
#define MM_SPARSE_STR  "coordinate"
#define MM_COMPLEX_STR  "complex"
#define MM_REAL_STR    "real"
#define MM_INT_STR    "integer"
#define MM_GENERAL_STR  "general"
#define MM_SYMM_STR    "symmetric"
#define MM_HERM_STR    "hermitian"
#define MM_SKEW_STR    "skew-symmetric"
#define MM_PATTERN_STR  "pattern"

#define mm_set_matrix(typecode)  ((*typecode)[0]='M')
#define mm_set_coordinate(typecode)  ((*typecode)[1]='C')
#define mm_set_array(typecode)  ((*typecode)[1]='A')
#define mm_set_dense(typecode)  mm_set_array(typecode)
#define mm_set_sparse(typecode)  mm_set_coordinate(typecode)

#define mm_set_complex(typecode)((*typecode)[2]='C')
#define mm_set_real(typecode)  ((*typecode)[2]='R')
#define mm_set_pattern(typecode)((*typecode)[2]='P')
#define mm_set_integer(typecode)((*typecode)[2]='I')

#define mm_set_symmetric(typecode)((*typecode)[3]='S')
#define mm_set_general(typecode)((*typecode)[3]='G')
#define mm_set_skew(typecode)  ((*typecode)[3]='K')
#define mm_set_hermitian(typecode)((*typecode)[3]='H')

#define mm_clear_typecode(typecode) ((*typecode)[0]=(*typecode)[1]= \
                  (*typecode)[2]=' ',(*typecode)[3]='G')

#define mm_initialize_typecode(typecode) mm_clear_typecode(typecode)

#define MM_COULD_NOT_READ_FILE  11
#define MM_PREMATURE_EOF    12
#define MM_NOT_MTX        13
#define MM_NO_HEADER      14
#define MM_UNSUPPORTED_TYPE    15
#define MM_LINE_TOO_LONG    16
#define MM_COULD_NOT_WRITE_FILE  17

typedef char MM_typecode[4];

static char *mm_typecode_to_str(MM_typecode matcode)
{
  char buffer[MM_MAX_LINE_LENGTH];
  char *types[4];
  int error =0;

  // check for MTX type
  if (mm_is_matrix(matcode))
    types[0] = MM_MTX_STR;
  else
    error=1;

  // check for CRD or ARR matrix
  if (mm_is_sparse(matcode))
    types[1] = MM_SPARSE_STR;
  else
    if (mm_is_dense(matcode))
      types[1] = MM_DENSE_STR;
    else
      return NULL;

  // check for element data type
  if (mm_is_real(matcode))
    types[2] = MM_REAL_STR;
  else
    if (mm_is_complex(matcode))
      types[2] = MM_COMPLEX_STR;
    else
      if (mm_is_pattern(matcode))
        types[2] = MM_PATTERN_STR;
      else
        if (mm_is_integer(matcode))
          types[2] = MM_INT_STR;
        else
          return NULL;

  // check for symmetry type
  if (mm_is_general(matcode))
    types[3] = MM_GENERAL_STR;
  else
    if (mm_is_symmetric(matcode))
      types[3] = MM_SYMM_STR;
    else
      if (mm_is_hermitian(matcode))
        types[3] = MM_HERM_STR;
      else
        if (mm_is_skew(matcode))
          types[3] = MM_SKEW_STR;
        else
          return NULL;

  sprintf(buffer,"%s %s %s %s", types[0], types[1], types[2], types[3]);
  return string_dup(buffer);
}

static int mm_read_banner(FILE *f, MM_typecode *matcode)
{
  char line[MM_MAX_LINE_LENGTH];
  char banner[MM_MAX_TOKEN_LENGTH];
  char mtx[MM_MAX_TOKEN_LENGTH];
  char crd[MM_MAX_TOKEN_LENGTH];
  char data_type[MM_MAX_TOKEN_LENGTH];
  char storage_scheme[MM_MAX_TOKEN_LENGTH];
  char *p;

  mm_clear_typecode(matcode);

  if (fgets(line, MM_MAX_LINE_LENGTH, f) == NULL)
    return MM_PREMATURE_EOF;

  if (sscanf(line, "%s %s %s %s %s", banner, mtx, crd, data_type,
        storage_scheme) != 5)
    return MM_PREMATURE_EOF;

  for (p=mtx; *p!='\0'; *p=(char)tolower(*p),p++); // convert to lower case
  for (p=crd; *p!='\0'; *p=(char)tolower(*p),p++);
  for (p=data_type; *p!='\0'; *p=(char)tolower(*p),p++);
  for (p=storage_scheme; *p!='\0'; *p=(char)tolower(*p),p++);

  // Check for banner.
  if (strncmp(banner, MatrixMarketBanner, strlen(MatrixMarketBanner)) != 0)
    return MM_NO_HEADER;

  // first field should be "mtx"
  if (strcmp(mtx, MM_MTX_STR) != 0)
    return  MM_UNSUPPORTED_TYPE;
  mm_set_matrix(matcode);

  // Second field describes whether this is a sparse matrix (in coordinate
  // storage) or a dense array .

  if (strcmp(crd, MM_SPARSE_STR) == 0)
    mm_set_sparse(matcode);
  else
    if (strcmp(crd, MM_DENSE_STR) == 0)
      mm_set_dense(matcode);
    else
      return MM_UNSUPPORTED_TYPE;


  // Third field.

  if (strcmp(data_type, MM_REAL_STR) == 0)
    mm_set_real(matcode);
  else
    if (strcmp(data_type, MM_COMPLEX_STR) == 0)
      mm_set_complex(matcode);
    else
      if (strcmp(data_type, MM_PATTERN_STR) == 0)
        mm_set_pattern(matcode);
      else
        if (strcmp(data_type, MM_INT_STR) == 0)
          mm_set_integer(matcode);
        else
          return MM_UNSUPPORTED_TYPE;


  // fourth field

  if (strcmp(storage_scheme, MM_GENERAL_STR) == 0)
    mm_set_general(matcode);
  else
  {
    if (strcmp(storage_scheme, MM_SYMM_STR) == 0)
      mm_set_symmetric(matcode);
    else
    {
      if (strcmp(storage_scheme, MM_HERM_STR) == 0)
        mm_set_hermitian(matcode);
      else
      {
        if (strcmp(storage_scheme, MM_SKEW_STR) == 0)
          mm_set_skew(matcode);
        else
          return MM_UNSUPPORTED_TYPE;
      }
    }
  }

  return 0;
}

static int mm_read_mtx_crd_size(FILE *f, size_t *M, size_t *N, size_t *nnz)
{
  char line[MM_MAX_LINE_LENGTH];
  int num_items_read;

  // set return null parameter values, in case we exit with errors
  *M = *N = *nnz = 0;

  // now continue scanning until you reach the end-of-comments
  do
  {
    if (fgets(line,MM_MAX_LINE_LENGTH,f) == NULL)
      return MM_PREMATURE_EOF;
  }
  while (line[0] == '%');

  // line[] is either blank or has M,N, nnz
  int iM, iN, innz;
  if (sscanf(line, "%d %d %d", &iM, &iN, &innz) == 3)
  {
    *M = (size_t)iM;
    *N = (size_t)iN;
    *nnz = (size_t)innz;
    return 0;
  }
  else
  {
    do
    {
      num_items_read = fscanf(f, "%d %d %d", &iM, &iN, &innz);
      if (num_items_read == EOF) return MM_PREMATURE_EOF;
    }
    while (num_items_read != 3);
    *M = (size_t)iM;
    *N = (size_t)iN;
    *nnz = (size_t)innz;
  }

  return 0;
}

static int mm_read_mtx_array_size(FILE *f, int *M, int *N)
{
  char line[MM_MAX_LINE_LENGTH];
  int num_items_read;

  // Set 0 parameter values, in case we exit with errors.
  *M = *N = 0;

  // Now continue scanning until you reach the end-of-comments.
  do
  {
    if (fgets(line,MM_MAX_LINE_LENGTH,f) == NULL)
      return MM_PREMATURE_EOF;
  }
  while (line[0] == '%');

  // line[] is either blank or has M,N, nz.
  if (sscanf(line, "%d %d", M, N) == 2)
    return 0;
  else // we have a blank line
  {
    do
    {
      num_items_read = fscanf(f, "%d %d", M, N);
      if (num_items_read == EOF) return MM_PREMATURE_EOF;
    }
    while (num_items_read != 2);
  }

  return 0;
}

static int mm_is_valid(MM_typecode matcode)
{
  if (!mm_is_matrix(matcode))
    return false;
  if (mm_is_dense(matcode) && mm_is_pattern(matcode))
    return false;
  if (mm_is_real(matcode) && mm_is_hermitian(matcode))
    return false;
  if (mm_is_pattern(matcode) &&
      (mm_is_hermitian(matcode) || mm_is_skew(matcode)))
    return false;
  return true;
}
//------------------------------------------------------------------------

// In-place distribution of global matrix.
static void distribute_matrix(krylov_factory_t* factory,
                              krylov_matrix_t** A,
                              matrix_sparsity_t** sparsity,
                              MPI_Comm comm,
                              index_t* row_dist)
{
  ASSERT(comm != MPI_COMM_SELF);

  int nprocs, rank;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);

  // Create a distributed sparsity pattern.
  distribute_matrix_sparsity(sparsity, comm, row_dist);

  // Create a distributed representation of the matrix.
  krylov_matrix_t* dist_A = krylov_factory_matrix(factory, *sparsity);

  // Shuttle the coefficients from the global into the local matrix.
  size_t nnz = matrix_sparsity_num_nonzeros(*sparsity);
  size_t num_local_rows = matrix_sparsity_num_local_rows(*sparsity);
  index_t rows[num_local_rows], cols[nnz];
  size_t num_cols[num_local_rows];
  int k = 0, rpos = 0, r = 0;
  index_t row;
  while (matrix_sparsity_next_row(*sparsity, &rpos, &row))
  {
    rows[r] = row;
    num_cols[r] = matrix_sparsity_num_columns(*sparsity, row);
    int pos = 0;
    index_t col;
    while (matrix_sparsity_next_column(*sparsity, row, &pos, &col))
      cols[k++] = col;
    ++r;
  }
  ASSERT(k == nnz);
  real_t values[nnz];
  krylov_matrix_get_values(*A, num_local_rows, num_cols, rows, cols, values);
  krylov_matrix_set_values(dist_A, num_local_rows, num_cols, rows, cols, values);

  // Set A to dist_A.
  krylov_matrix_free(*A);
  *A = dist_A;
}

static krylov_matrix_t* krylov_factory_matrix_from_mm(krylov_factory_t* factory,
                                                      MPI_Comm comm,
                                                      FILE* f)
{
  // Rewind the file descriptor in case we've used it.
  fseek(f, 0, SEEK_SET);

  // Read the banner to get the type code.
  MM_typecode matcode;
  mm_read_banner(f, &matcode);

  // We only support sparse matrices.
  if (!mm_is_sparse(matcode))
  {
    polymec_error("krylov_factory_matrix_from_mm: Not a sparse matrix type: %s",
                  mm_typecode_to_str(matcode));
  }

  // We don't support complex coefficients.
  if (mm_is_complex(matcode))
  {
    polymec_error("krylov_factory_matrix_from_mm: unsupported matrix type: %s",
                  mm_typecode_to_str(matcode));
  }

  // Size the thing up.
  size_t M, N, nnz;
  if (mm_read_mtx_crd_size(f, &M, &N, &nnz) != 0)
    polymec_error("krylov_factory_matrix_from_mm: error reading matrix size.");

  if (M != N)
    polymec_error("krylov_factory_matrix_from_mm: only square matrices are supported.");

  // Dream up a naive partitioning.
  int nprocs;
  MPI_Comm_size(comm, &nprocs);
  index_t row_dist[nprocs+1];
  row_dist[0] = 0;
  for (int p = 0; p < nprocs-1; ++p)
    row_dist[p+1] = row_dist[p] + (int)(1.0*M/nprocs);
  row_dist[nprocs] = M;

  // Retrieve the values into arrays.
  index_t* I = polymec_malloc(nnz * sizeof(index_t));
  index_t* J = polymec_malloc(nnz * sizeof(index_t));
  real_t* val = polymec_malloc(nnz * sizeof(real_t));
  index_t num_cols[M];
  memset(num_cols, 0, M * sizeof(index_t));
  for (index_t i = 0; i < nnz; ++i)
  {
#if POLYMEC_HAVE_SINGLE_PRECISION
    int num_items = fscanf(f, "%" SCNu64 " %" SCNu64 " %g\n", &I[i], &J[i], &val[i]);
#else
    int num_items = fscanf(f, "%" SCNu64 " %" SCNu64 " %lg\n", &I[i], &J[i], &val[i]);
#endif
    if (num_items != 3)
      goto error;

    // adjust from 1-based to 0-based
    --I[i];
    --J[i];

    // Tally the columns.
    ++num_cols[I[i]];
  }

  // Construct a global sparsity pattern.
  index_t global_row_dist[2] = {0, M};
  matrix_sparsity_t* sparsity = matrix_sparsity_new(MPI_COMM_SELF, global_row_dist);
  {
    // Set up the number of columns in each row.
    int rpos = 0, r = 0;
    index_t row;
    while (matrix_sparsity_next_row(sparsity, &rpos, &row))
    {
      matrix_sparsity_set_num_columns(sparsity, row, num_cols[r]);
      ++r;
    }

    // Set up column counters for placing edges in the graph.
    index_t* columns[M], column_counter[M];
    for (int i = 0; i < M; ++i)
    {
      columns[i] = matrix_sparsity_columns(sparsity, i);
      column_counter[i] = 0;
    }

    // Stick all the columns in.
    for (index_t k = 0; k < nnz; ++k)
    {
      index_t Ik = I[k], Jk = J[k];
      columns[Ik][column_counter[Ik]] = Jk;
      ++column_counter[Ik];
    }
  }

  // Create a fresh matrix.
  krylov_matrix_t* A = krylov_factory_matrix(factory, sparsity);

  // Insert the values into the matrix. This is dumb and slow, but simple.
  // And we've already hit the disk, so it's not a big deal.
  for (int i = 0; i < nnz; ++i)
  {
    size_t num_columns = 1;
    index_t row = I[i];
    index_t column = J[i];
    krylov_matrix_set_values(A, 1, &num_columns, &row, &column, &val[i]);
  }
  krylov_matrix_assemble(A);

  // Clean up.
  polymec_free(I);
  polymec_free(J);
  polymec_free(val);

  // Distribute as necessary.
  if ((comm != MPI_COMM_SELF) && (nprocs > 1))
  {
    distribute_matrix(factory, &A, &sparsity, comm, row_dist);
    matrix_sparsity_free(sparsity);
    return A;
  }
  else
  {
    matrix_sparsity_free(sparsity);
    return A;
  }

error:
  polymec_free(I);
  polymec_free(J);
  polymec_free(val);
  polymec_error("krylov_factory_matrix_from_mm: invalid data read.");
  return NULL;
}

krylov_matrix_t* krylov_factory_matrix_from_file(krylov_factory_t* factory,
                                                 MPI_Comm comm,
                                                 const char* filename)
{
  FILE* f = fopen(filename, "r");

  // Make sure the file exists.
  if (f == NULL)
    polymec_error("krylov_factory_matrix_from_file: Could not read file %s.", filename);

  krylov_matrix_t* A = NULL;

  // Make a guess that the format of the matrix and dispatch as necessary.
  MM_typecode matcode;
  if ((mm_read_banner(f, &matcode) == 0) && (mm_is_valid(matcode)))
    A = krylov_factory_matrix_from_mm(factory, comm, f);
  else
    polymec_error("krylov_factory_matrix_from_file: unsupported format");

  fclose(f);

  krylov_matrix_assemble(A);

  return A;
}

void krylov_matrix_free(krylov_matrix_t* A)
{
  if ((A->context != NULL) && (A->vtable.dtor != NULL))
    A->vtable.dtor(A->context);
  polymec_free(A);
}

MPI_Comm krylov_matrix_comm(krylov_matrix_t* A)
{
  return A->comm;
}

size_t krylov_matrix_block_size(krylov_matrix_t* A, index_t block_row)
{
  if (A->vtable.block_size != NULL)
    return A->vtable.block_size(A->context, block_row);
  else
    return 1;
}

krylov_matrix_t* krylov_matrix_clone(krylov_matrix_t* A)
{
  krylov_matrix_t* B = polymec_malloc(sizeof(krylov_matrix_t));
  B->context = A->vtable.clone(A->context);
  B->vtable = A->vtable;
  B->comm = A->comm;
  B->num_local_rows = A->num_local_rows;
  B->num_global_rows = A->num_global_rows;
  return B;
}

void krylov_matrix_copy(krylov_matrix_t* A, krylov_matrix_t* copy)
{
  START_FUNCTION_TIMER();
  ASSERT(copy->num_global_rows == A->num_global_rows);
  ASSERT(copy->num_local_rows == A->num_local_rows);
  A->vtable.copy(A->context, copy->context);
  STOP_FUNCTION_TIMER();
}

void* krylov_matrix_impl(krylov_matrix_t* A)
{
  return A->context;
}

size_t krylov_matrix_num_local_rows(krylov_matrix_t* A)
{
  return A->num_local_rows;
}

size_t krylov_matrix_num_global_rows(krylov_matrix_t* A)
{
  return A->num_global_rows;
}

void krylov_matrix_zero(krylov_matrix_t* A)
{
  START_FUNCTION_TIMER();
  A->vtable.zero(A->context);
  STOP_FUNCTION_TIMER();
}

void krylov_matrix_add_identity(krylov_matrix_t* A,
                                real_t scale_factor)
{
  START_FUNCTION_TIMER();
  A->vtable.add_identity(A->context, scale_factor);
  STOP_FUNCTION_TIMER();
}

void krylov_matrix_scale(krylov_matrix_t* A,
                         real_t scale_factor)
{
  START_FUNCTION_TIMER();
  A->vtable.scale(A->context, scale_factor);
  STOP_FUNCTION_TIMER();
}

void krylov_matrix_diag_scale(krylov_matrix_t* A,
                              krylov_vector_t* L,
                              krylov_vector_t* R)
{
  START_FUNCTION_TIMER();
  A->vtable.diag_scale(A->context, L->context, R->context);
  STOP_FUNCTION_TIMER();
}

void krylov_matrix_add_diagonal(krylov_matrix_t* A,
                                krylov_vector_t* D)
{
  A->vtable.add_diagonal(A->context, D->context);
}

void krylov_matrix_set_diagonal(krylov_matrix_t* A,
                                krylov_vector_t* D)
{
  START_FUNCTION_TIMER();
  A->vtable.set_diagonal(A->context, D->context);
  STOP_FUNCTION_TIMER();
}

void krylov_matrix_matvec(krylov_matrix_t* A,
                          krylov_vector_t* x,
                          bool transpose,
                          krylov_vector_t* y)
{
  START_FUNCTION_TIMER();
  A->vtable.matvec(A->context, x->context, transpose, y->context);
  STOP_FUNCTION_TIMER();
}

void krylov_matrix_set_values(krylov_matrix_t* A,
                              size_t num_rows,
                              size_t* num_columns,
                              index_t* rows, index_t* columns,
                              real_t* values)
{
  START_FUNCTION_TIMER();
  A->vtable.set_values(A->context, num_rows, num_columns, rows, columns, values);
  STOP_FUNCTION_TIMER();
}

void krylov_matrix_add_values(krylov_matrix_t* A,
                              size_t num_rows,
                              size_t* num_columns,
                              index_t* rows, index_t* columns,
                              real_t* values)
{
  START_FUNCTION_TIMER();
  A->vtable.add_values(A->context, num_rows, num_columns, rows, columns, values);
  STOP_FUNCTION_TIMER();
}

void krylov_matrix_set_blocks(krylov_matrix_t* A,
                              size_t num_blocks,
                              index_t* block_rows,
                              index_t* block_columns,
                              real_t* block_values)
{
  START_FUNCTION_TIMER();
  if (A->vtable.set_blocks != NULL)
    A->vtable.set_blocks(A->context, num_blocks, block_rows, block_columns, block_values);
  else
    polymec_error("Non-block matrix cannot use block interface.");
  STOP_FUNCTION_TIMER();
}

void krylov_matrix_add_blocks(krylov_matrix_t* A,
                              size_t num_blocks,
                              index_t* block_rows,
                              index_t* block_columns,
                              real_t* block_values)
{
  START_FUNCTION_TIMER();
  if (A->vtable.add_blocks != NULL)
    A->vtable.add_blocks(A->context, num_blocks, block_rows, block_columns, block_values);
  else
    polymec_error("Non-block matrix cannot use block interface.");
  STOP_FUNCTION_TIMER();
}

void krylov_matrix_get_blocks(krylov_matrix_t* A,
                              size_t num_blocks,
                              index_t* block_rows,
                              index_t* block_columns,
                              real_t* block_values)
{
  START_FUNCTION_TIMER();
  if (A->vtable.add_blocks != NULL)
    A->vtable.get_blocks(A->context, num_blocks, block_rows, block_columns, block_values);
  else
    polymec_error("Non-block matrix cannot use block interface.");
  STOP_FUNCTION_TIMER();
}

void krylov_matrix_set_block(krylov_matrix_t* A,
                             index_t block_row,
                             index_t block_column,
                             real_t* block_values)
{
  START_FUNCTION_TIMER();
  krylov_matrix_set_blocks(A, 1, &block_row, &block_column, block_values);
  STOP_FUNCTION_TIMER();
}

void krylov_matrix_add_block(krylov_matrix_t* A,
                             index_t block_row,
                             index_t block_column,
                             real_t* block_values)
{
  krylov_matrix_add_blocks(A, 1, &block_row, &block_column, block_values);
}

void krylov_matrix_get_block(krylov_matrix_t* A,
                             index_t block_row,
                             index_t block_column,
                             real_t* block_values)
{
  krylov_matrix_get_blocks(A, 1, &block_row, &block_column, block_values);
}

void krylov_matrix_assemble(krylov_matrix_t* A)
{
  START_FUNCTION_TIMER();
  if (A->vtable.assemble != NULL)
    A->vtable.assemble(A->context);
  STOP_FUNCTION_TIMER();
}

void krylov_matrix_get_values(krylov_matrix_t* A,
                              size_t num_rows,
                              size_t* num_columns,
                              index_t* rows,
                              index_t* columns,
                              real_t* values)
{
  START_FUNCTION_TIMER();
  A->vtable.get_values(A->context, num_rows, num_columns, rows, columns, values);
  STOP_FUNCTION_TIMER();
}

void krylov_matrix_fprintf(krylov_matrix_t* A,
                           FILE* stream)
{
  if ((A->vtable.fprintf != NULL) && (stream != NULL))
    A->vtable.fprintf(A->context, stream);
}

//------------------------------------------------------------------------
//                          Krylov vector
//------------------------------------------------------------------------

krylov_vector_t* krylov_vector_new(void* context,
                                   krylov_vector_vtable vtable,
                                   size_t local_size,
                                   size_t global_size)
{
  ASSERT(vtable.clone != NULL);
  ASSERT(vtable.copy != NULL);
  ASSERT(vtable.zero != NULL);
  ASSERT(vtable.set_value != NULL);
  ASSERT(vtable.scale != NULL);
  ASSERT(vtable.diag_scale != NULL);
  ASSERT(vtable.set_values != NULL);
  ASSERT(vtable.add_values != NULL);
  ASSERT(vtable.get_values != NULL);
  ASSERT(vtable.copy_in != NULL);
  ASSERT(vtable.copy_out != NULL);
  ASSERT(vtable.dot != NULL);
  ASSERT(vtable.norm != NULL);
  ASSERT(vtable.w2_norm != NULL);
  ASSERT(vtable.wrms_norm != NULL);
  ASSERT(local_size > 0);
  ASSERT(global_size >= local_size);
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

void krylov_vector_copy(krylov_vector_t* v, krylov_vector_t* copy)
{
  START_FUNCTION_TIMER();
  ASSERT(copy->local_size == v->local_size);
  ASSERT(copy->global_size == v->global_size);
  v->vtable.copy(v->context, copy->context);
  STOP_FUNCTION_TIMER();
}

void* krylov_vector_impl(krylov_vector_t* v)
{
  return v->context;
}

size_t krylov_vector_local_size(krylov_vector_t* v)
{
  return v->local_size;
}

size_t krylov_vector_global_size(krylov_vector_t* v)
{
  return v->global_size;
}

void krylov_vector_zero(krylov_vector_t* v)
{
  START_FUNCTION_TIMER();
  v->vtable.zero(v->context);
  STOP_FUNCTION_TIMER();
}

void krylov_vector_set_value(krylov_vector_t* v,
                             real_t value)
{
  START_FUNCTION_TIMER();
  v->vtable.set_value(v->context, value);
  STOP_FUNCTION_TIMER();
}

void krylov_vector_scale(krylov_vector_t* v,
                         real_t scale_factor)
{
  START_FUNCTION_TIMER();
  v->vtable.scale(v->context, scale_factor);
  STOP_FUNCTION_TIMER();
}

void krylov_vector_diag_scale(krylov_vector_t* v,
                              krylov_vector_t* D)
{
  START_FUNCTION_TIMER();
  ASSERT(D != NULL);
  v->vtable.diag_scale(v->context, D->context);
  STOP_FUNCTION_TIMER();
}

void krylov_vector_set_values(krylov_vector_t* v,
                              size_t num_values,
                              index_t* indices,
                              real_t* values)
{
  START_FUNCTION_TIMER();
  v->vtable.set_values(v->context, num_values, indices, values);
  STOP_FUNCTION_TIMER();
}

void krylov_vector_add_values(krylov_vector_t* v,
                              size_t num_values,
                              index_t* indices,
                              real_t* values)
{
  START_FUNCTION_TIMER();
  v->vtable.add_values(v->context, num_values, indices, values);
  STOP_FUNCTION_TIMER();
}

void krylov_vector_get_values(krylov_vector_t* v,
                              size_t num_values,
                              index_t* indices,
                              real_t* values)
{
  START_FUNCTION_TIMER();
  v->vtable.get_values(v->context, num_values, indices, values);
  STOP_FUNCTION_TIMER();
}

void krylov_vector_copy_in(krylov_vector_t* v,
                           real_t* local_values)
{
  START_FUNCTION_TIMER();
  v->vtable.copy_in(v->context, local_values);
  STOP_FUNCTION_TIMER();
}

void krylov_vector_copy_out(krylov_vector_t* v,
                            real_t* local_values)
{
  START_FUNCTION_TIMER();
  v->vtable.copy_out(v->context, local_values);
  STOP_FUNCTION_TIMER();
}

void krylov_vector_assemble(krylov_vector_t* v)
{
  START_FUNCTION_TIMER();
  if (v->vtable.assemble != NULL)
    v->vtable.assemble(v->context);
  STOP_FUNCTION_TIMER();
}

real_t krylov_vector_dot(krylov_vector_t* v,
                         krylov_vector_t* w)
{
  START_FUNCTION_TIMER();
  ASSERT(v != NULL);
  ASSERT(w != NULL);
  real_t product = v->vtable.dot(v->context, w->context);
  STOP_FUNCTION_TIMER();
  return product;
}

real_t krylov_vector_norm(krylov_vector_t* v, int p)
{
  START_FUNCTION_TIMER();
  ASSERT((p == 0) || (p == 1) || (p == 2));
  real_t norm = v->vtable.norm(v->context, p);
  STOP_FUNCTION_TIMER();
  return norm;
}

real_t krylov_vector_w2_norm(krylov_vector_t* v, krylov_vector_t* w)
{
  START_FUNCTION_TIMER();
  real_t norm;
  if (w != NULL)
    norm = v->vtable.w2_norm(v->context, w->context);
  else
    norm = v->vtable.norm(v->context, 2);
  STOP_FUNCTION_TIMER();
  return norm;
}

real_t krylov_vector_wrms_norm(krylov_vector_t* v, krylov_vector_t* w)
{
  START_FUNCTION_TIMER();
  ASSERT(w != NULL);
  real_t norm = v->vtable.wrms_norm(v->context, w->context);
  STOP_FUNCTION_TIMER();
  return norm;
}

void krylov_vector_fprintf(krylov_vector_t* v,
                           FILE* stream)
{
  if ((v->vtable.fprintf != NULL) && (stream != NULL))
    v->vtable.fprintf(v->context, stream);
}

static void distribute_vector(krylov_factory_t* factory,
                              krylov_vector_t** x,
                              MPI_Comm comm,
                              index_t* row_dist)
{
  START_FUNCTION_TIMER();
  ASSERT(comm != MPI_COMM_SELF);
  ASSERT((*x)->local_size == (*x)->global_size);

  int nprocs, rank;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);

  // Create a local representation of the vector.
  krylov_vector_t* dist_x = krylov_factory_vector(factory, comm, row_dist);
  size_t num_local_rows = (size_t)(row_dist[rank+1] - row_dist[rank]);

  // Shuttle the values from the global into the local vector.
  index_t rows[num_local_rows];
  real_t values[num_local_rows];
  for (int i = 0; i < (int)num_local_rows; ++i)
  {
    index_t I = row_dist[rank] + i; // global vertex index
    rows[i] = I;
  }
  krylov_vector_get_values(*x, num_local_rows, rows, values);
  krylov_vector_set_values(dist_x, num_local_rows, rows, values);

  krylov_vector_free(*x);
  *x = dist_x;
  STOP_FUNCTION_TIMER();
}

static krylov_vector_t* krylov_factory_vector_from_mm(krylov_factory_t* factory,
                                                      MPI_Comm comm,
                                                      FILE* f)
{
  START_FUNCTION_TIMER();

  // Rewind the file descriptor in case we've used it.
  fseek(f, 0, SEEK_SET);

  // Read the banner to get the type code.
  MM_typecode matcode;
  mm_read_banner(f, &matcode);

  // A vector is an Nx1 dense matrix.
  if (!mm_is_dense(matcode))
  {
    polymec_error("krylov_factory_vector_from_mm: Not a vector/dense matrix type: %s",
                  mm_typecode_to_str(matcode));
  }

  // We don't support complex coefficients.
  if (mm_is_complex(matcode))
  {
    polymec_error("krylov_factory_vector_from_mm: unsupported vector type: %s",
                  mm_typecode_to_str(matcode));
  }

  // Size the thing up.
  int M, N;
  if (mm_read_mtx_array_size(f, &M, &N) != 0)
    polymec_error("krylov_factory_vector_from_mm: error reading vector size.");
  if (N != 1)
    polymec_error("krylov_factory_vector_from_mm: cannot create vector from %dx%d matrix.");

  // Dream up a naive partitioning.
  int nprocs;
  MPI_Comm_size(comm, &nprocs);
  index_t row_dist[nprocs+1];
  row_dist[0] = 0;
  for (int p = 0; p < nprocs-1; ++p)
    row_dist[p+1] = row_dist[p] + (int)(1.0*M/nprocs);
  row_dist[nprocs] = M;

  // Retrieve the values into an array.
  index_t rows[M];
  real_t values[M];
  size_t count = 0;
  for (int i = 0; i < M; ++i)
  {
    rows[i] = (index_t)i;
#if POLYMEC_HAVE_SINGLE_PRECISION
    count += fscanf(f, "%g\n", &values[i]);
#else
    count += fscanf(f, "%lg\n", &values[i]);
#endif
  }
  ASSERT(count == (size_t)M);

  index_t self_row_dist[2] = {0, M};
  krylov_vector_t* x = krylov_factory_vector(factory, MPI_COMM_SELF, self_row_dist);

  // Insert the values into the vector.
  krylov_vector_set_values(x, M, rows, values);
  krylov_vector_assemble(x);

  // Distribute as necessary.
  if (comm != MPI_COMM_SELF)
    distribute_vector(factory, &x, comm, row_dist);

  STOP_FUNCTION_TIMER();
  return x;
}

krylov_vector_t* krylov_factory_vector_from_file(krylov_factory_t* factory,
                                                 MPI_Comm comm,
                                                 const char* filename)
{
  START_FUNCTION_TIMER();

  // Make sure the file exists.
  FILE* f = fopen(filename, "r");
  if (f == NULL)
    polymec_error("krylov_factory_vector_from_file: Could not read file %s.", filename);

  krylov_vector_t* x = NULL;

  // Make a guess that the format of the vector and dispatch as necessary.
  MM_typecode matcode;
  if ((mm_read_banner(f, &matcode) == 0) && (mm_is_valid(matcode)))
    x = krylov_factory_vector_from_mm(factory, comm, f);
  else
    polymec_error("krylov_factory_vector_from_file: unsupported format");

  fclose(f);

  krylov_vector_assemble(x);

  STOP_FUNCTION_TIMER();
  return x;
}

//------------------------------------------------------------------------
//                  Factories for creating Krylov solvers
//------------------------------------------------------------------------

krylov_factory_t* krylov_factory_new(const char* name,
                                     void* context,
                                     krylov_factory_vtable vtable)
{
  ASSERT(vtable.pcg_solver != NULL);
  ASSERT(vtable.gmres_solver != NULL);
  ASSERT(vtable.bicgstab_solver != NULL);
  ASSERT(vtable.preconditioner != NULL);
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
                                       matrix_sparsity_t* sparsity)
{
  return factory->vtable.matrix(factory->context, sparsity);
}

krylov_matrix_t* krylov_factory_block_matrix(krylov_factory_t* factory,
                                             matrix_sparsity_t* sparsity,
                                             size_t block_size)
{
  ASSERT(block_size > 0);
  return factory->vtable.block_matrix(factory->context, sparsity, block_size);
}

krylov_matrix_t* krylov_factory_var_block_matrix(krylov_factory_t* factory,
                                                 matrix_sparsity_t* sparsity,
                                                 size_t* block_sizes)
{
  ASSERT(block_sizes != NULL);
  if (factory->vtable.var_block_matrix != NULL)
    return factory->vtable.var_block_matrix(factory->context, sparsity, block_sizes);
  else
    polymec_fatal_error("The %s factory does not support variable block matrices.", factory->name);
}

krylov_vector_t* krylov_factory_vector(krylov_factory_t* factory,
                                       MPI_Comm comm,
                                       index_t* row_distribution)
{
  return factory->vtable.vector(factory->context, comm, row_distribution);
}

krylov_solver_t* krylov_factory_pcg_solver(krylov_factory_t* factory,
                                           MPI_Comm comm)
{
  return factory->vtable.pcg_solver(factory->context, comm);
}

krylov_solver_t* krylov_factory_gmres_solver(krylov_factory_t* factory,
                                             MPI_Comm comm,
                                             int krylov_dimension)
{
  ASSERT(krylov_dimension >= 3);
  return factory->vtable.gmres_solver(factory->context, comm, krylov_dimension);
}

krylov_solver_t* krylov_factory_bicgstab_solver(krylov_factory_t* factory,
                                                MPI_Comm comm)
{
  return factory->vtable.bicgstab_solver(factory->context, comm);
}

krylov_solver_t* krylov_factory_special_solver(krylov_factory_t* factory,
                                               MPI_Comm comm,
                                               const char* solver_name,
                                               string_string_unordered_map_t* options)
{
  if (factory->vtable.special_solver != NULL)
    return factory->vtable.special_solver(factory->context, comm, solver_name, options);
  else
    return NULL;
}

krylov_pc_t* krylov_factory_preconditioner(krylov_factory_t* factory,
                                           MPI_Comm comm,
                                           const char* pc_name,
                                           string_string_unordered_map_t* options)
{
  return factory->vtable.preconditioner(factory->context, comm, pc_name, options);
}

