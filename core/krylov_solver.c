// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <float.h>
#include "core/krylov_solver.h"
#include "core/timer.h"
#include "core/string_utils.h"

//------------------------------------------------------------------------
//                Krylov data types and virtual tables
//------------------------------------------------------------------------

struct krylov_solver_t
{
  char* name;
  void* context;
  krylov_solver_vtable vtable;
  krylov_pc_t* pc;
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
  int num_local_rows;
  int num_global_rows;
};

struct krylov_vector_t
{
  void* context;
  krylov_vector_vtable vtable;
  int local_size;
  index_t global_size;
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
  krylov_solver_t* solver = polymec_malloc(sizeof(krylov_solver_t));
  solver->name = string_dup(name);
  solver->context = context;
  solver->vtable = vtable;
  solver->pc = NULL;
  return solver;
}

void krylov_solver_free(krylov_solver_t* solver)
{
  if ((solver->context != NULL) && (solver->vtable.dtor != NULL))
    solver->vtable.dtor(solver->context);
  if (solver->pc != NULL)
    krylov_pc_free(solver->pc);
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
  return solver->vtable.solve(solver->context, b->context, x->context, 
                              residual_norm, num_iterations);
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
                                   int num_local_rows,
                                   index_t num_global_rows)
{
  ASSERT(vtable.zero != NULL);
  ASSERT(vtable.add_diagonal != NULL);
  ASSERT(vtable.set_diagonal != NULL);
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
// This code was pulled from http://math.nist.gov/MatrixMarket.
//------------------------------------------------------------------------

#define MM_MAX_LINE_LENGTH 1025
#define MatrixMarketBanner "%%MatrixMarket"
#define MM_MAX_TOKEN_LENGTH 64

#define mm_is_matrix(typecode)	((typecode)[0]=='M')

#define mm_is_sparse(typecode)	((typecode)[1]=='C')
#define mm_is_coordinate(typecode)((typecode)[1]=='C')
#define mm_is_dense(typecode)	((typecode)[1]=='A')
#define mm_is_array(typecode)	((typecode)[1]=='A')

#define mm_is_complex(typecode)	((typecode)[2]=='C')
#define mm_is_real(typecode)		((typecode)[2]=='R')
#define mm_is_pattern(typecode)	((typecode)[2]=='P')
#define mm_is_integer(typecode) ((typecode)[2]=='I')

#define mm_is_symmetric(typecode)((typecode)[3]=='S')
#define mm_is_general(typecode)	((typecode)[3]=='G')
#define mm_is_skew(typecode)	((typecode)[3]=='K')
#define mm_is_hermitian(typecode)((typecode)[3]=='H')

#define MM_MTX_STR		"matrix"
#define MM_ARRAY_STR	"array"
#define MM_DENSE_STR	"array"
#define MM_COORDINATE_STR "coordinate" 
#define MM_SPARSE_STR	"coordinate"
#define MM_COMPLEX_STR	"complex"
#define MM_REAL_STR		"real"
#define MM_INT_STR		"integer"
#define MM_GENERAL_STR  "general"
#define MM_SYMM_STR		"symmetric"
#define MM_HERM_STR		"hermitian"
#define MM_SKEW_STR		"skew-symmetric"
#define MM_PATTERN_STR  "pattern"

#define mm_set_matrix(typecode)	((*typecode)[0]='M')
#define mm_set_coordinate(typecode)	((*typecode)[1]='C')
#define mm_set_array(typecode)	((*typecode)[1]='A')
#define mm_set_dense(typecode)	mm_set_array(typecode)
#define mm_set_sparse(typecode)	mm_set_coordinate(typecode)

#define mm_set_complex(typecode)((*typecode)[2]='C')
#define mm_set_real(typecode)	((*typecode)[2]='R')
#define mm_set_pattern(typecode)((*typecode)[2]='P')
#define mm_set_integer(typecode)((*typecode)[2]='I')

#define mm_set_symmetric(typecode)((*typecode)[3]='S')
#define mm_set_general(typecode)((*typecode)[3]='G')
#define mm_set_skew(typecode)	((*typecode)[3]='K')
#define mm_set_hermitian(typecode)((*typecode)[3]='H')

#define mm_clear_typecode(typecode) ((*typecode)[0]=(*typecode)[1]= \
									(*typecode)[2]=' ',(*typecode)[3]='G')

#define mm_initialize_typecode(typecode) mm_clear_typecode(typecode)

#define MM_COULD_NOT_READ_FILE	11
#define MM_PREMATURE_EOF		12
#define MM_NOT_MTX				13
#define MM_NO_HEADER			14
#define MM_UNSUPPORTED_TYPE		15
#define MM_LINE_TOO_LONG		16
#define MM_COULD_NOT_WRITE_FILE	17

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

  for (p=mtx; *p!='\0'; *p=tolower(*p),p++); // convert to lower case 
  for (p=crd; *p!='\0'; *p=tolower(*p),p++);  
  for (p=data_type; *p!='\0'; *p=tolower(*p),p++);
  for (p=storage_scheme; *p!='\0'; *p=tolower(*p),p++);

  // check for banner 
  if (strncmp(banner, MatrixMarketBanner, strlen(MatrixMarketBanner)) != 0)
    return MM_NO_HEADER;

  // first field should be "mtx" 
  if (strcmp(mtx, MM_MTX_STR) != 0)
    return  MM_UNSUPPORTED_TYPE;
  mm_set_matrix(matcode);

  // second field describes whether this is a sparse matrix (in coordinate
  // storage) or a dense array 

  if (strcmp(crd, MM_SPARSE_STR) == 0)
    mm_set_sparse(matcode);
  else
    if (strcmp(crd, MM_DENSE_STR) == 0)
      mm_set_dense(matcode);
    else
      return MM_UNSUPPORTED_TYPE;


  // third field 

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

static int mm_read_mtx_crd_size(FILE *f, int *M, int *N, int *nz)
{
  char line[MM_MAX_LINE_LENGTH];
  int num_items_read;

  // set return null parameter values, in case we exit with errors 
  *M = *N = *nz = 0;

  // now continue scanning until you reach the end-of-comments 
  do 
  {
    if (fgets(line,MM_MAX_LINE_LENGTH,f) == NULL) 
      return MM_PREMATURE_EOF;
  }
  while (line[0] == '%');

  // line[] is either blank or has M,N, nz 
  if (sscanf(line, "%d %d %d", M, N, nz) == 3)
    return 0;

  else
  {
    do
    { 
      num_items_read = fscanf(f, "%d %d %d", M, N, nz); 
      if (num_items_read == EOF) return MM_PREMATURE_EOF;
    }
    while (num_items_read != 3);
  }

  return 0;
}

static int mm_read_mtx_array_size(FILE *f, int *M, int *N)
{
  char line[MM_MAX_LINE_LENGTH];
  int num_items_read;
  // set return null parameter values, in case we exit with errors 
  *M = *N = 0;

  // now continue scanning until you reach the end-of-comments 
  do 
  {
    if (fgets(line,MM_MAX_LINE_LENGTH,f) == NULL) 
      return MM_PREMATURE_EOF;
  }
  while (line[0] == '%');

  // line[] is either blank or has M,N, nz 
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
  if (!mm_is_matrix(matcode)) return false;
  if (mm_is_dense(matcode) && mm_is_pattern(matcode)) return false;
  if (mm_is_real(matcode) && mm_is_hermitian(matcode)) return false;
  if (mm_is_pattern(matcode) && (mm_is_hermitian(matcode) || 
        mm_is_skew(matcode))) return false;
  return true;
}

static int mm_read_mtx_crd_data(FILE *f, int M, int N, int nz, 
                                int I[], int J[], real_t val[], 
                                MM_typecode matcode)
{
  int i;
  if (mm_is_complex(matcode))
  {
    for (i=0; i<nz; i++)
      if (fscanf(f, "%d %d %lg %lg", &I[i], &J[i], &val[2*i], &val[2*i+1])
          != 4) return MM_PREMATURE_EOF;
  }
  else if (mm_is_real(matcode))
  {
    for (i=0; i<nz; i++)
    {
      if (fscanf(f, "%d %d %lg\n", &I[i], &J[i], &val[i])
          != 3) return MM_PREMATURE_EOF;

    }
  }

  else if (mm_is_pattern(matcode))
  {
    for (i=0; i<nz; i++)
      if (fscanf(f, "%d %d", &I[i], &J[i])
          != 2) return MM_PREMATURE_EOF;
  }
  else
    return MM_UNSUPPORTED_TYPE;

  return 0;
}

static int mm_read_unsymmetric_sparse(const char *fname, 
                                      int *M_, int *N_, int *nz_,
                                      real_t **val_, int **I_, int **J_)
{

  FILE *f = fopen(fname, "r");
  if (f == NULL)
    return -1;

  MM_typecode matcode;
  if (mm_read_banner(f, &matcode) != 0)
  {
    printf("mm_read_unsymmetric: Could not process Matrix Market banner ");
    printf(" in file [%s]\n", fname);
    return -1;
  }

  if (!(mm_is_real(matcode) && mm_is_matrix(matcode) &&
       mm_is_sparse(matcode)))
  {
    fprintf(stderr, "Sorry, this application does not support ");
    fprintf(stderr, "Market Market type: [%s]\n",
        mm_typecode_to_str(matcode));
    return -1;
  }

  // find out size of sparse matrix: M, N, nz .... 
  int M, N, nz;
  if (mm_read_mtx_crd_size(f, &M, &N, &nz) !=0)
  {
    fprintf(stderr, "read_unsymmetric_sparse(): could not parse matrix size.\n");
    return -1;
  }

  *M_ = M;
  *N_ = N;
  *nz_ = nz;

  // reserve memory for matrices 
  int* I = polymec_malloc(nz * sizeof(int));
  int* J = polymec_malloc(nz * sizeof(int));
  real_t* val = polymec_malloc(nz * sizeof(real_t));

  *val_ = val;
  *I_ = I;
  *J_ = J;

  // NOTE: when reading in reals, ANSI C requires the use of the "l"  
  //   specifier as in "%lg", "%lf", "%le", otherwise errors will occur 
  //  (ANSI C X3.159-1989, Sec. 4.9.6.2, p. 136 lines 13-15)            

  for (int i = 0; i < nz; ++i)
  {
#if POLYMEC_HAVE_SINGLE_PRECISION
    fscanf(f, "%d %d %g\n", &I[i], &J[i], &val[i]);
#else
    fscanf(f, "%d %d %lg\n", &I[i], &J[i], &val[i]);
#endif
    I[i]--;  // adjust from 1-based to 0-based 
    J[i]--;
  }
  fclose(f);

  return 0;
}
//------------------------------------------------------------------------
static krylov_matrix_t* redistribute_matrix(krylov_factory_t* factory, 
                                            krylov_matrix_t* global_A, 
                                            adj_graph_t* global_sparsity, 
                                            MPI_Comm comm)
{
  ASSERT(comm != MPI_COMM_SELF);

  int nprocs, rank;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_size(comm, &rank);

  // Do a naive partitioning (since we assume a sensible partitioning was 
  // done to create the row space for the matrix in the first place...)
  int num_global_rows = adj_graph_num_vertices(global_sparsity);
  int num_local_rows = num_global_rows / nprocs;
  if (rank == (nprocs - 1))
    num_local_rows = num_global_rows - rank * num_local_rows;

  // Break up the graph.
  adj_graph_t* local_sparsity = adj_graph_new(comm, num_local_rows);
  int tot_num_edges = 0;
  for (int i = 0; i < num_local_rows; ++i)
  {
    int I = rank * num_local_rows + i; // global vertex index
    int num_edges = adj_graph_num_edges(global_sparsity, I);
    int* global_edges = adj_graph_edges(global_sparsity, I);
    adj_graph_set_num_edges(local_sparsity, i, num_edges);
    int* local_edges = adj_graph_edges(local_sparsity, i);
    for (int j = 0; j < num_edges; ++j)
      local_edges[j] = global_edges[j];
    tot_num_edges += num_edges;
  }

  // Create a local representation of the matrix.
  krylov_matrix_t* local_A = krylov_factory_matrix(factory, local_sparsity);
  adj_graph_free(local_sparsity);

  // Shuttle the coefficients from the global into the local matrix.
  index_t global_rows[num_local_rows], local_rows[num_local_rows], 
          num_cols[num_local_rows], cols[tot_num_edges];
  real_t global_vals[tot_num_edges];
  int k = 0;
  for (int i = 0; i < num_local_rows; ++i)
  {
    local_rows[i] = i;
    index_t I = rank * num_local_rows + i; // global vertex index
    global_rows[i] = I;
    num_cols[i] = adj_graph_num_edges(global_sparsity, I);
    int* global_edges = adj_graph_edges(global_sparsity, I);
    for (int j = 0; j < num_cols[i]; ++j, ++k)
      cols[++k] = global_edges[j];
  }
  ASSERT(k == tot_num_edges);
  real_t values[tot_num_edges];
  krylov_matrix_get_values(global_A, num_local_rows, num_cols, global_rows, 
                           cols, values);
  krylov_matrix_set_values(local_A, num_local_rows, num_cols, local_rows, 
                           cols, values);
  return local_A;
}

static krylov_matrix_t* krylov_factory_matrix_from_mm(krylov_factory_t* factory,
                                                      MPI_Comm comm,
                                                      const char* filename)
{
  FILE* f = fopen(filename, "r");
  ASSERT(f != NULL);

  // Read the banner to get the type code.
  MM_typecode matcode;
  if (mm_read_banner(f, &matcode) != 0)
    polymec_error("krylov_factory_matrix_from_mm: %s is not in Matrix Market format.", filename);

  // We don't support complex coefficients.
  if (mm_is_complex(matcode))
  {
    polymec_error("krylov_factory_matrix_from_mm: unsupported matrix type: %s", 
                  mm_typecode_to_str(matcode));
  }

  // Size the thing up.
  int M, N, nz;
  if (mm_read_mtx_crd_size(f, &M, &N, &nz) != 0)
    polymec_error("krylov_factory_matrix_from_mm: error reading matrix size.");

  // Retrieve the values into a packed data structure.
  int* I = polymec_malloc(nz * sizeof(int));
  int* J = polymec_malloc(nz * sizeof(int));
  real_t* val = polymec_malloc(nz * sizeof(real_t));
  int num_offdiags[M];
  memset(num_offdiags, 0, M * sizeof(int));
  for (int i = 0; i < nz; ++i)
  {
#if POLYMEC_HAVE_SINGLE_PRECISION
    fscanf(f, "%d %d %g\n", &I[i], &J[i], &val[i]);
#else
    fscanf(f, "%d %d %lg\n", &I[i], &J[i], &val[i]);
#endif
    // adjust from 1-based to 0-based 
    I[i]--;  
    J[i]--;

    // Tally the off-diagonal elements.
    if (J[i] != I[i])
      ++num_offdiags[I[i]];
  }

  // We're through with the file.
  fclose(f);

  // Construct a local sparsity graph.
  adj_graph_t* sparsity = adj_graph_new(MPI_COMM_WORLD, M);
  for (int i = 0; i < M; ++i)
    adj_graph_set_num_edges(sparsity, i, num_offdiags[i]);
  int* edges[M], edge_counter[M];
  for (int i = 0; i < M; ++i)
  {
    edges[i] = adj_graph_edges(sparsity, i);
    edge_counter[i] = 0;
  }
  for (int i = 0; i < nz; ++i)
  {
    edges[i][edge_counter[i]] = J[i];
    ++edge_counter[i];
  }

  // Create a fresh matrix.
  krylov_matrix_t* global_A = krylov_factory_matrix(factory, sparsity);

  // Insert the values into the matrix. This is dumb and slow, but simple.
  // And we've already hit the disk, so it's not a big deal.
  for (int i = 0; i < nz; ++i)
  {
    index_t num_columns = 1;
    index_t row = I[i];
    index_t column = J[i];
    krylov_matrix_set_values(global_A, 1, &num_columns, &row, &column, &val[i]);
  }

  // Clean up.
  polymec_free(I);
  polymec_free(J);
  polymec_free(val);

  // Distribute as necessary.
  if (comm != MPI_COMM_SELF)
  {
    krylov_matrix_t* local_A = redistribute_matrix(factory, global_A, sparsity, comm);
    krylov_matrix_free(global_A);
    adj_graph_free(sparsity);
    return local_A;
  }
  else
  {
    adj_graph_free(sparsity);
    return global_A;
  }
}

krylov_matrix_t* krylov_factory_matrix_from_file(krylov_factory_t* factory,  
                                                 MPI_Comm comm,
                                                 const char* filename,
                                                 krylov_matrix_format_t format)
{
  // Make sure the file exists.
  if (!file_exists(filename))
    polymec_error("krylov_factory_matrix_from_file: Could not read file %s.", filename);

  krylov_matrix_t* global_A;
  if (format == MATRIX_MARKET)
  {
    return krylov_factory_matrix_from_mm(factory, comm, filename);
  }
  else
    polymec_error("krylov_factory_matrix_from_file: unsupported format");
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
                                krylov_vector_t* D)
{
  A->vtable.add_diagonal(A->context, D->context);
}

void krylov_matrix_set_diagonal(krylov_matrix_t* A,
                                krylov_vector_t* D)
{
  A->vtable.set_diagonal(A->context, D->context);
}

void krylov_matrix_set_values(krylov_matrix_t* A,
                              index_t num_rows,
                              index_t* num_columns,
                              index_t* rows, index_t* columns,
                              real_t* values)
{
  A->vtable.set_values(A->context, num_rows, num_columns, rows, columns, values);
}
                              
void krylov_matrix_add_values(krylov_matrix_t* A,
                              index_t num_rows,
                              index_t* num_columns,
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
  if (A->vtable.start_assembly != NULL)
    A->vtable.start_assembly(A->context);
}

void krylov_matrix_finish_assembly(krylov_matrix_t* A)
{
  if (A->vtable.finish_assembly != NULL)
    A->vtable.finish_assembly(A->context);
}

void krylov_matrix_get_values(krylov_matrix_t* A,
                              index_t num_rows,
                              index_t* num_columns,
                              index_t* rows, index_t* columns,
                              real_t* values)
{
  A->vtable.get_values(A->context, num_rows, num_columns, rows, columns, values);
}

//------------------------------------------------------------------------
//                          Krylov vector
//------------------------------------------------------------------------

krylov_vector_t* krylov_vector_new(void* context,
                                   krylov_vector_vtable vtable,
                                   int local_size,
                                   index_t global_size)
{
  ASSERT(vtable.zero != NULL);
  ASSERT(vtable.set_value != NULL);
  ASSERT(vtable.scale != NULL);
  ASSERT(vtable.set_values != NULL);
  ASSERT(vtable.add_values != NULL);
  ASSERT(vtable.get_values != NULL);
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
                              index_t num_values,
                              index_t* indices,
                              real_t* values)
{
  v->vtable.set_values(v->context, num_values, indices, values);
}
                              
void krylov_vector_add_values(krylov_vector_t* v,
                              index_t num_values,
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
  if (v->vtable.start_assembly != NULL)
    v->vtable.start_assembly(v->context);
}

void krylov_vector_finish_assembly(krylov_vector_t* v)
{
  if (v->vtable.finish_assembly != NULL)
    v->vtable.finish_assembly(v->context);
}

void krylov_vector_get_values(krylov_vector_t* v,
                              index_t num_values,
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

krylov_matrix_t* krylov_factory_var_block_matrix(krylov_factory_t* factory, 
                                                 adj_graph_t* sparsity,
                                                 int* block_sizes)
{
  ASSERT(block_sizes != NULL);
  if (factory->vtable.var_block_matrix != NULL)
    return factory->vtable.var_block_matrix(factory->context, sparsity, block_sizes);
  else
    polymec_error("The %s factory does not support variable block matrices.", factory->name);
}

krylov_vector_t* krylov_factory_vector(krylov_factory_t* factory,
                                       adj_graph_t* dist_graph)
{
  return factory->vtable.vector(factory->context, dist_graph);
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

