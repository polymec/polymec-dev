// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>

#include "cmocka.h"
#include "core/polymec.h"
#include "core/file_utils.h"
#include "core/krylov_solver.h"

// This helper creates a sparsity pattern for a 1D finite difference 
// discretization of the Laplacian operator.
static matrix_sparsity_t* create_1d_laplacian_sparsity(MPI_Comm comm, int N_local)
{
  int nprocs, rank;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);

  index_t row_dist[nprocs+1];
  row_dist[0] = 0;
  for (int p = 0; p < nprocs; ++p)
    row_dist[p+1] = row_dist[p] + N_local;
  matrix_sparsity_t* sparsity = matrix_sparsity_new(comm, row_dist);
  index_t N_global = matrix_sparsity_num_global_rows(sparsity);

  int rpos = 0;
  index_t row;
  while (matrix_sparsity_next_row(sparsity, &rpos, &row))
  {
    if (row == 0)
    {
      matrix_sparsity_set_num_columns(sparsity, row, 2);
      index_t* cols = matrix_sparsity_columns(sparsity, row);
      cols[0] = row;
      cols[1] = row+1;
    }
    else if (row == (N_global-1))
    {
      matrix_sparsity_set_num_columns(sparsity, row, 2);
      index_t* cols = matrix_sparsity_columns(sparsity, row);
      cols[0] = row-1;
      cols[1] = row;
    }
    else
    {
      matrix_sparsity_set_num_columns(sparsity, row, 3);
      index_t* cols = matrix_sparsity_columns(sparsity, row);
      cols[0] = row-1;
      cols[1] = row;
      cols[2] = row+1;
    }
  }

  return sparsity;
}

// This helper creates a sparsity pattern for a 2D finite difference 
// discretization of the Laplacian operator.
static matrix_sparsity_t* create_2d_laplacian_sparsity(MPI_Comm comm, 
                                                       int Nx_local,
                                                       int Ny)
{
  int nprocs, rank;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);

  index_t row_dist[nprocs+1];
  row_dist[0] = 0;
  for (int p = 0; p < nprocs; ++p)
    row_dist[p+1] = row_dist[p] + Ny * Nx_local;
  matrix_sparsity_t* sparsity = matrix_sparsity_new(comm, row_dist);
  index_t N_global = matrix_sparsity_num_global_rows(sparsity);

  int rpos = 0;
  index_t row;
  while (matrix_sparsity_next_row(sparsity, &rpos, &row))
  {
    if (row == 0)
    {
      matrix_sparsity_set_num_columns(sparsity, row, 2*Ny);
      index_t* cols = matrix_sparsity_columns(sparsity, row);
      for (int c = 0; c < 2*Ny; ++c)
        cols[c] = row+c;
    }
    else if (row == (N_global-1))
    {
      matrix_sparsity_set_num_columns(sparsity, row, 2*Ny);
      index_t* cols = matrix_sparsity_columns(sparsity, row);
      for (int c = 0; c < 2*Ny; ++c)
        cols[c] = row-2*Ny+c+1;
    }
    else
    {
      matrix_sparsity_set_num_columns(sparsity, row, 3*Ny);
      index_t* cols = matrix_sparsity_columns(sparsity, row);
      for (int c = 0; c < 3*Ny; ++c)
        cols[c] = row-2*Ny+c+1;
    }
  }

  return sparsity;
}

#if POLYMEC_HAVE_SHARED_LIBS
static krylov_factory_t* create_petsc_krylov_factory()
{
  krylov_factory_t* factory = NULL;
  char* petsc_dir = getenv("PETSC_DIR");
  if (petsc_dir != NULL)
  {
    char* petsc_arch = getenv("PETSC_ARCH");
    if (petsc_arch != NULL)
    {
      // Check for the existence of the PETSc dynamic library.
      char petsc_path[FILENAME_MAX+1];
      snprintf(petsc_path, FILENAME_MAX, "%s/%s/lib/libpetsc%s", petsc_dir, petsc_arch, SHARED_LIBRARY_SUFFIX);
      if (file_exists(petsc_path))
      {
        factory = petsc_krylov_factory(petsc_dir, petsc_arch);
        if (factory == NULL)
          log_urgent("Could not load PETSc. Skipping PETSc test.");
      }
      else
        log_urgent("PETSc library not found. Skipping PETSc test.");
    }
    else
      log_urgent("PETSC_ARCH not set. Skipping PETSc test.");
  }
  else
    log_urgent("PETSC_DIR not set. Skipping PETSc test.");
  return factory;
}

static krylov_factory_t* create_hypre_krylov_factory()
{
  krylov_factory_t* factory = NULL;
  char* hypre_dir = getenv("HYPRE_DIR"); // Absent other info, we rely on this.
  if (hypre_dir != NULL)
  {
    // Check for the existence of the HYPRE dynamic library.
    char hypre_path[FILENAME_MAX+1];
    snprintf(hypre_path, FILENAME_MAX, "%s/libHYPRE%s", hypre_dir, SHARED_LIBRARY_SUFFIX);
    if (file_exists(hypre_path))
    {
      factory = hypre_krylov_factory(hypre_dir);
      if (factory == NULL)
        log_urgent("Could not load HYPRE. Skipping HYPRE test.");
    }
    else
      log_urgent("HYPRE library not found. Skipping HYPRE test.");
  }
  else
    log_urgent("HYPRE_DIR not set. Skipping HYPRE test.");
  return factory;
}
#endif

static void test_krylov_factory(void** state, krylov_factory_t* factory)
{
  if (factory != NULL)
  {
    assert_true(krylov_factory_name(factory) != NULL);
    krylov_factory_free(factory);
  }
}

static void test_krylov_matrix(void** state, krylov_factory_t* factory)
{
  if (factory != NULL)
  {
    MPI_Comm comm = MPI_COMM_WORLD;

    // Create a sparsity pattern with 1000 local rows.
    matrix_sparsity_t* sparsity = create_1d_laplacian_sparsity(comm, 1000);

    // Create a matrix.
    krylov_matrix_t* mat = krylov_factory_matrix(factory, sparsity);
    assert_int_equal(1000, krylov_matrix_num_local_rows(mat));
    int nprocs;
    MPI_Comm_size(comm, &nprocs);
    assert_int_equal(1000*nprocs, krylov_matrix_num_global_rows(mat));

    // Clone it.
    krylov_matrix_t* mat1 = krylov_matrix_clone(mat);
    assert_int_equal(1000, krylov_matrix_num_local_rows(mat1));

    // Put everything away.
    krylov_matrix_free(mat);
    krylov_matrix_free(mat1);
    krylov_factory_free(factory);
    matrix_sparsity_free(sparsity);
  }
}

static void test_krylov_matrix_from_file(void** state, 
                                         krylov_factory_t* factory, 
                                         const char* filename,
                                         int num_global_rows,
                                         int num_test_values,
                                         index_t* test_rows,
                                         index_t* test_cols,
                                         real_t* test_values)
{
  if (factory != NULL)
  {
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank, nprocs;
    MPI_Comm_size(comm, &nprocs);
    MPI_Comm_rank(comm, &rank);
    int num_local_rows = num_global_rows / nprocs;
    if (rank == (nprocs-1))
      num_local_rows = num_global_rows - rank*num_local_rows;

    // Read a matrix in from the given file and evaluate its data.
    krylov_matrix_t* mat = krylov_factory_matrix_from_file(factory, comm, filename);
    assert_int_equal(num_local_rows, krylov_matrix_num_local_rows(mat));
    assert_int_equal(num_global_rows, krylov_matrix_num_global_rows(mat));
    for (int i = 0; i < num_test_values; ++i)
    {
      // Only test locally-maintained values.
      index_t r = test_rows[i] - rank*num_local_rows;
      if (r < num_local_rows)
      {
        index_t one = 1;
        real_t val;
        krylov_matrix_get_values(mat, 1, &one, 
                                 &test_rows[i], &test_cols[i], &val);
        assert_true(ABS(val - test_values[i]) < 1e-14);
      }
    }
    krylov_matrix_free(mat);
    krylov_factory_free(factory);
  }
}

static void test_krylov_vector(void** state, krylov_factory_t* factory)
{
  if (factory != NULL)
  {
    MPI_Comm comm = MPI_COMM_WORLD;
    int nprocs;
    MPI_Comm_size(comm, &nprocs);

    // Create a vector with 1000 rows per process.
    index_t N_local = 1000;
    index_t row_dist[nprocs+1];
    row_dist[0] = 0;
    for (int p = 0; p < nprocs; ++p)
      row_dist[p+1] = row_dist[p] + N_local;
    krylov_vector_t* vec = krylov_factory_vector(factory, comm, row_dist);
    assert_int_equal(1000, krylov_vector_local_size(vec));
    assert_int_equal(1000*nprocs, krylov_vector_global_size(vec));

    // Clone it.
    krylov_vector_t* vec1 = krylov_vector_clone(vec);
    assert_int_equal(1000, krylov_vector_local_size(vec1));

    // Put everything away.
    krylov_vector_free(vec);
    krylov_vector_free(vec1);
    krylov_factory_free(factory);
  }
}

static void test_krylov_vector_from_file(void** state, 
                                         krylov_factory_t* factory, 
                                         const char* filename,
                                         int num_global_rows,
                                         int num_test_values,
                                         index_t* test_indices,
                                         real_t* test_values)
{
  if (factory != NULL)
  {
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank, nprocs;
    MPI_Comm_size(comm, &nprocs);
    MPI_Comm_rank(comm, &rank);
    int num_local_rows = num_global_rows / nprocs;
    if (rank == (nprocs-1))
      num_local_rows = num_global_rows - rank*num_local_rows;

    // Read a vector in from the given file.
    krylov_vector_t* vec = krylov_factory_vector_from_file(factory, comm, filename);
    assert_true(krylov_vector_local_size(vec) > 0);
    assert_int_equal(num_local_rows, krylov_vector_local_size(vec));
    assert_int_equal(num_global_rows, krylov_vector_global_size(vec));
    for (int i = 0; i < num_test_values; ++i)
    {
      // Only test locally-maintained values.
      index_t r = test_indices[i] - rank*num_local_rows;
      if (r < num_local_rows)
      {
        real_t val;
        krylov_vector_get_values(vec, 1, &test_indices[i], &val);
        assert_true(ABS(val - test_values[i]) < 1e-14);
      }
    }
    krylov_vector_free(vec);
    krylov_factory_free(factory);
  }
}

static void test_krylov_matrix_from_sherman1(void** state, krylov_factory_t* factory)
{
  int num_test_values = 4;
  index_t test_rows[] = {0, 14, 22, 238};
  index_t test_cols[] = {0, 4, 22, 228};
  real_t test_values[] = {-0.005649, 1.127e-5, -2.541, 0.1127};
  test_krylov_matrix_from_file(state, 
                               factory, 
                               CMAKE_CURRENT_SOURCE_DIR "/sherman1.mtx",
                               1000,
                               num_test_values, test_rows, test_cols, test_values);
}

static void test_krylov_vector_from_sherman1_b(void** state, krylov_factory_t* factory)
{
  int num_test_values = 4;
  index_t test_indices[] = {4, 294, 295, 491};
  real_t test_values[] = {2.254e-5, -9.148e-10, 0.0, -2.429e-9};
  test_krylov_vector_from_file(state, 
                               factory, 
                               CMAKE_CURRENT_SOURCE_DIR "/sherman1_b.mtx",
                               1000,
                               num_test_values, test_indices, test_values);
}

typedef enum
{
  PCG_SOLVER,
  GMRES_SOLVER,
  BICGSTAB_SOLVER
} iterative_solver_t;

static void test_1d_laplace_eqn(void** state, 
                                krylov_factory_t* factory,
                                iterative_solver_t solver_type)
{
  if (factory != NULL)
  {
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank, nprocs;
    MPI_Comm_size(comm, &nprocs);
    MPI_Comm_rank(comm, &rank);

    // Create a distributed matrix with 100 rows per process.
    int N = 100;
    matrix_sparsity_t* sparsity = create_1d_laplacian_sparsity(comm, N);
    index_t* row_dist = matrix_sparsity_row_distribution(sparsity);
    index_t N_global = matrix_sparsity_num_global_rows(sparsity);
    krylov_matrix_t* A = krylov_factory_matrix(factory, sparsity);

    // Create a RHS vector.
    krylov_vector_t* b = krylov_factory_vector(factory, comm, row_dist);

    // Construct the 1D Laplacian operator and RHS.
    real_t h = 1.0 / N;
    int rpos = 0;
    index_t row;
    while (matrix_sparsity_next_row(sparsity, &rpos, &row))
    {
      real_t Aij[3], bi[1];
      index_t rows[3], cols[3], num_cols;
      if (row == 0)
      {
        num_cols = 2; 
        rows[0] = 0, cols[0] = 0, Aij[0] = -2.0;
                     cols[1] = 1, Aij[1] = 1.0;
        krylov_matrix_set_values(A, 1, &num_cols, rows, cols, Aij);
        bi[0] = 0.0;
        krylov_vector_set_values(b, 1, &row, bi);
      }
      else if (row == (N_global-1))
      {
        num_cols = 2;
        rows[0] = row, cols[0] = row-1, Aij[0] = 1.0;
                       cols[1] = row,   Aij[1] = -2.0;
        krylov_matrix_set_values(A, 1, &num_cols, rows, cols, Aij);
        bi[0] = 1.0;
        krylov_vector_set_values(b, 1, &row, bi);
      }
      else
      {
        num_cols = 3;
        rows[0] = row, cols[0] = row-1, Aij[0] = 1.0;
        rows[1] = row, cols[1] = row,   Aij[1] = -2.0;
        rows[2] = row, cols[2] = row+1, Aij[2] = 1.0;
        krylov_matrix_set_values(A, 1, &num_cols, rows, cols, Aij);
        bi[0] = 0.0;
        krylov_vector_set_values(b, 1, &row, bi);
      }
    }
    krylov_matrix_assemble(A);
    krylov_vector_assemble(b);
    krylov_matrix_scale(A, 1.0/(h*h));
    krylov_vector_scale(b, -1.0/(h*h));
//krylov_matrix_fprintf(A, stdout);
//krylov_vector_fprintf(b, stdout);

    // Create a solution vector.
    krylov_vector_t* x = krylov_factory_vector(factory, comm, row_dist);

    // Create a solver.
    krylov_solver_t* solver;
    switch (solver_type)
    {
      case PCG_SOLVER: solver = krylov_factory_pcg_solver(factory, comm); break;
      case GMRES_SOLVER: solver = krylov_factory_gmres_solver(factory, comm, 30); break;
      case BICGSTAB_SOLVER: solver = krylov_factory_bicgstab_solver(factory, comm); break;
    }
    assert_true(solver != NULL);
    krylov_solver_set_tolerances(solver, 1e-5, 1e-8, 2.0);
    if (solver_type == GMRES_SOLVER)
      krylov_solver_set_max_iterations(solver, nprocs * 1000);
    else
      krylov_solver_set_max_iterations(solver, 1000);

    // Solve the equation.
    krylov_solver_set_operator(solver, A);
    real_t res_norm;
    int num_iters;
    bool solved = krylov_solver_solve(solver, b, x, &res_norm, &num_iters);
    log_debug("residual norm is %g, # iterations is %d", res_norm, num_iters);
    assert_true(solved);

    if (log_level() == LOG_DEBUG)
    {
      printf("x =\n");
      krylov_vector_fprintf(x, stdout);
    }

    // Put everything away.
    krylov_solver_free(solver);
    krylov_matrix_free(A);
    krylov_vector_free(x);
    krylov_vector_free(b);
    krylov_factory_free(factory);
  }
}

static void test_load_and_solve(void** state, 
                                krylov_factory_t* factory, 
                                const char* mat_filename,
                                const char* rhs_filename)
{
  if (factory != NULL)
  {
    MPI_Comm comm = MPI_COMM_WORLD;

    // Read in a matrix A from the given file.
    krylov_matrix_t* A = krylov_factory_matrix_from_file(factory, comm, mat_filename);

    // Read in a vector b from the given file.
    krylov_vector_t* b = krylov_factory_vector_from_file(factory, comm, rhs_filename);

    // Create a solution vector.
    krylov_vector_t* x = krylov_vector_clone(b);
    krylov_vector_zero(x);

    // Create a solver.
    krylov_solver_t* solver = krylov_factory_bicgstab_solver(factory, comm);
    assert_true(solver != NULL);
    krylov_solver_set_tolerances(solver, 1e-5, 1e-8, 2.0);
    krylov_solver_set_max_iterations(solver, 100);

    // Solve the equation.
    krylov_solver_set_operator(solver, A);
    real_t res_norm;
    int num_iters;
    bool solved = krylov_solver_solve(solver, b, x, &res_norm, &num_iters);
    log_debug("residual norm is %g, # iterations is %d", res_norm, num_iters);
    assert_true(solved);

    // Clean up.
    krylov_solver_free(solver);
    krylov_matrix_free(A);
    krylov_vector_free(x);
    krylov_vector_free(b);
    krylov_factory_free(factory);
  }
}

static void test_2d_laplace_eqn(void** state, 
                                krylov_factory_t* factory,
                                iterative_solver_t solver_type,
                                bool use_block_matrix)
{
  if (factory != NULL)
  {
    MPI_Comm comm = MPI_COMM_WORLD;
    int rank, nprocs;
    MPI_Comm_size(comm, &nprocs);
    MPI_Comm_rank(comm, &rank);

    // Create a distributed matrix with 10 block rows per process, 
    // with 10x10 blocks T.
    int N = 10;
    static real_t T[10*10] = {-4.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                               1.0,-4.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                               0.0, 1.0,-4.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                               0.0, 0.0, 1.0,-4.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                               0.0, 0.0, 0.0, 1.0,-4.0, 1.0, 0.0, 0.0, 0.0, 0.0, 
                               0.0, 0.0, 0.0, 0.0, 1.0,-4.0, 1.0, 0.0, 0.0, 0.0, 
                               0.0, 0.0, 0.0, 0.0, 0.0, 1.0,-4.0, 1.0, 0.0, 0.0, 
                               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,-4.0, 1.0, 0.0, 
                               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,-4.0, 1.0, 
                               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,-4.0};
    static real_t b_left[10] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    static real_t b_mid[10] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0};
    static real_t b_right[10] = {-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0,-1.0};
    matrix_sparsity_t* sparsity;
    if (use_block_matrix)
      sparsity = create_1d_laplacian_sparsity(comm, N);
    else
      sparsity = create_2d_laplacian_sparsity(comm, N*N, N*N);
    index_t* row_dist = matrix_sparsity_row_distribution(sparsity);
    index_t N_global = matrix_sparsity_num_global_rows(sparsity);
    krylov_matrix_t* A;
    if (use_block_matrix)
      A = krylov_factory_block_matrix(factory, sparsity, N);
    else
      A = krylov_factory_matrix(factory, sparsity);

    // Create a RHS vector.
    krylov_vector_t* b = krylov_factory_vector(factory, comm, row_dist);

    // Construct the 2D Laplacian operator and RHS.
    real_t h = 1.0 / N;
    int rpos = 0;
    index_t row;
    if (use_block_matrix)
    {
      // The sparsity pattern uses block rows, so we insert blocks 
      // for each row.
      while (matrix_sparsity_next_row(sparsity, &rpos, &row))
      {
        if (row == 0)
        {
          krylov_matrix_set_block(A, 0, 0, T);
          krylov_matrix_set_block(A, 0, 1, T);
          index_t rows[N];
          for (int i = 0; i < N; ++i)
            rows[i] = i;
          krylov_vector_set_values(b, N, rows, b_left);
        }
        else if (row == (N_global-1))
        {
          krylov_matrix_set_block(A, N_global-1, N_global-2, T);
          krylov_matrix_set_block(A, N_global-1, N_global-1, T);
          index_t rows[N];
          for (int i = 0; i < N; ++i)
            rows[i] = N_global - N + i;
          krylov_vector_set_values(b, N, rows, b_right);
        }
        else
        {
          krylov_matrix_set_block(A, row, row-1, T);
          krylov_matrix_set_block(A, row, row, T);
          krylov_matrix_set_block(A, row, row+1, T);
          index_t rows[N];
          for (int i = 0; i < N; ++i)
            rows[i] = N*row + i;
          krylov_vector_set_values(b, N, rows, b_mid);
        }
      }
    }
    else
    {
      // Our matrix sparsity doesn't use block rows, but we still want to 
      // insert blocks into the matrix to show that this works, so we 
      // loop over block rows.
      for (index_t row = row_dist[rank]/N; row < row_dist[rank+1]/N; ++row)
      {
        if (row == 0)
        {
          krylov_matrix_set_block(A, 0, 0, T);
          krylov_matrix_set_block(A, 0, 1, T);
          index_t rows[N];
          for (int i = 0; i < N; ++i)
            rows[i] = i;
          krylov_vector_set_values(b, N, rows, b_left);
        }
        else if (row == (N_global-1))
        {
          krylov_matrix_set_block(A, N_global-1, N_global-2, T);
          krylov_matrix_set_block(A, N_global-1, N_global-1, T);
          index_t rows[N];
          for (int i = 0; i < N; ++i)
            rows[i] = N_global - N + i;
          krylov_vector_set_values(b, N, rows, b_right);
        }
        else
        {
          krylov_matrix_set_block(A, row, row-1, T);
          krylov_matrix_set_block(A, row, row, T);
          krylov_matrix_set_block(A, row, row+1, T);
          index_t rows[N];
          for (int i = 0; i < N; ++i)
            rows[i] = N*row + i;
          krylov_vector_set_values(b, N, rows, b_mid);
        }
      }
    }
    krylov_matrix_assemble(A);
    krylov_vector_assemble(b);
    krylov_matrix_scale(A, 1.0/(h*h));
    krylov_vector_scale(b, -1.0/(h*h));
//krylov_matrix_fprintf(A, stdout);
//krylov_vector_fprintf(b, stdout);

    // Create a solution vector.
    krylov_vector_t* x = krylov_factory_vector(factory, comm, row_dist);

    // Create a solver.
    krylov_solver_t* solver;
    switch (solver_type)
    {
      case PCG_SOLVER: solver = krylov_factory_pcg_solver(factory, comm); break;
      case GMRES_SOLVER: solver = krylov_factory_gmres_solver(factory, comm, 30); break;
      case BICGSTAB_SOLVER: solver = krylov_factory_bicgstab_solver(factory, comm); break;
    }
    assert_true(solver != NULL);
    krylov_solver_set_tolerances(solver, 1e-5, 1e-8, 2.0);
    if (solver_type == GMRES_SOLVER)
      krylov_solver_set_max_iterations(solver, nprocs * 1000);
    else
      krylov_solver_set_max_iterations(solver, 1000);

    // Solve the equation.
    krylov_solver_set_operator(solver, A);
    real_t res_norm;
    int num_iters;
    bool solved = krylov_solver_solve(solver, b, x, &res_norm, &num_iters);
    log_debug("residual norm is %g, # iterations is %d", res_norm, num_iters);
    assert_true(solved);

    if (log_level() == LOG_DEBUG)
    {
      printf("x =\n");
      krylov_vector_fprintf(x, stdout);
    }

    // Put everything away.
    krylov_solver_free(solver);
    krylov_matrix_free(A);
    krylov_vector_free(x);
    krylov_vector_free(b);
    krylov_factory_free(factory);
  }
}

void test_lis_krylov_factory(void** state)
{
  krylov_factory_t* lis = lis_krylov_factory();
  test_krylov_factory(state, lis);
}

void test_lis_krylov_matrix(void** state)
{
  krylov_factory_t* lis = lis_krylov_factory();
  test_krylov_matrix(state, lis);
}

void test_lis_krylov_matrix_from_file(void** state)
{
  krylov_factory_t* lis = lis_krylov_factory();
  test_krylov_matrix_from_sherman1(state, lis);
}

void test_lis_krylov_vector(void** state)
{
  krylov_factory_t* lis = lis_krylov_factory();
  test_krylov_vector(state, lis);
}

void test_lis_krylov_vector_from_file(void** state)
{
  krylov_factory_t* lis = lis_krylov_factory();
  test_krylov_vector_from_sherman1_b(state, lis);
}

void test_lis_pcg_1d_laplace_eqn(void** state)
{
  krylov_factory_t* lis = lis_krylov_factory();
  test_1d_laplace_eqn(state, lis, PCG_SOLVER);
}

void test_lis_gmres_1d_laplace_eqn(void** state)
{
  krylov_factory_t* lis = lis_krylov_factory();
  test_1d_laplace_eqn(state, lis, GMRES_SOLVER);
}

void test_lis_bicgstab_1d_laplace_eqn(void** state)
{
  krylov_factory_t* lis = lis_krylov_factory();
  test_1d_laplace_eqn(state, lis, BICGSTAB_SOLVER);
}

void test_lis_sherman1(void** state)
{
  krylov_factory_t* lis = lis_krylov_factory();
  test_load_and_solve(state, lis, 
                      CMAKE_CURRENT_SOURCE_DIR "/sherman1.mtx", 
                      CMAKE_CURRENT_SOURCE_DIR "/sherman1_b.mtx");
}

void test_lis_pcg_2d_laplace_eqn(void** state)
{
  krylov_factory_t* lis = lis_krylov_factory();
  test_2d_laplace_eqn(state, lis, PCG_SOLVER, false);
}

void test_lis_pcg_block_2d_laplace_eqn(void** state)
{
  krylov_factory_t* lis = lis_krylov_factory();
  test_2d_laplace_eqn(state, lis, PCG_SOLVER, true);
}

void test_lis_gmres_2d_laplace_eqn(void** state)
{
  krylov_factory_t* lis = lis_krylov_factory();
  test_2d_laplace_eqn(state, lis, GMRES_SOLVER, false);
}

void test_lis_gmres_block_2d_laplace_eqn(void** state)
{
  krylov_factory_t* lis = lis_krylov_factory();
  test_2d_laplace_eqn(state, lis, GMRES_SOLVER, true);
}

void test_lis_bicgstab_2d_laplace_eqn(void** state)
{
  krylov_factory_t* lis = lis_krylov_factory();
  test_2d_laplace_eqn(state, lis, BICGSTAB_SOLVER, false);
}

void test_lis_bicgstab_block_2d_laplace_eqn(void** state)
{
  krylov_factory_t* lis = lis_krylov_factory();
  test_2d_laplace_eqn(state, lis, BICGSTAB_SOLVER, true);
}

#if POLYMEC_HAVE_SHARED_LIBS
void test_petsc_krylov_factory(void** state)
{
  krylov_factory_t* petsc = create_petsc_krylov_factory();
  test_krylov_factory(state, petsc);
}

void test_petsc_krylov_matrix(void** state)
{
  krylov_factory_t* petsc = create_petsc_krylov_factory();
  test_krylov_matrix(state, petsc);
}

void test_petsc_krylov_matrix_from_file(void** state)
{
  krylov_factory_t* petsc = create_petsc_krylov_factory();
  test_krylov_matrix_from_sherman1(state, petsc);
}

void test_petsc_krylov_vector(void** state)
{
  krylov_factory_t* petsc = create_petsc_krylov_factory();
  test_krylov_vector(state, petsc);
}

void test_petsc_krylov_vector_from_file(void** state)
{
  krylov_factory_t* petsc = create_petsc_krylov_factory();
  test_krylov_vector_from_sherman1_b(state, petsc);
}

void test_petsc_pcg_1d_laplace_eqn(void** state)
{
  krylov_factory_t* petsc = create_petsc_krylov_factory();
  test_1d_laplace_eqn(state, petsc, PCG_SOLVER);
}

void test_petsc_gmres_1d_laplace_eqn(void** state)
{
  krylov_factory_t* petsc = create_petsc_krylov_factory();
  test_1d_laplace_eqn(state, petsc, GMRES_SOLVER);
}

void test_petsc_bicgstab_1d_laplace_eqn(void** state)
{
  krylov_factory_t* petsc = create_petsc_krylov_factory();
  test_1d_laplace_eqn(state, petsc, BICGSTAB_SOLVER);
}

void test_petsc_sherman1(void** state)
{
  krylov_factory_t* petsc = create_petsc_krylov_factory();
  test_load_and_solve(state, petsc, 
                      CMAKE_CURRENT_SOURCE_DIR "/sherman1.mtx", 
                      CMAKE_CURRENT_SOURCE_DIR "/sherman1_b.mtx");
}

void test_petsc_pcg_2d_laplace_eqn(void** state)
{
  krylov_factory_t* petsc = create_petsc_krylov_factory();
  test_2d_laplace_eqn(state, petsc, PCG_SOLVER, false);
}

void test_petsc_pcg_block_2d_laplace_eqn(void** state)
{
  krylov_factory_t* petsc = create_petsc_krylov_factory();
  test_2d_laplace_eqn(state, petsc, PCG_SOLVER, true);
}

void test_petsc_gmres_2d_laplace_eqn(void** state)
{
  krylov_factory_t* petsc = create_petsc_krylov_factory();
  test_2d_laplace_eqn(state, petsc, GMRES_SOLVER, false);
}

void test_petsc_gmres_block_2d_laplace_eqn(void** state)
{
  krylov_factory_t* petsc = create_petsc_krylov_factory();
  test_2d_laplace_eqn(state, petsc, GMRES_SOLVER, true);
}

void test_petsc_bicgstab_2d_laplace_eqn(void** state)
{
  krylov_factory_t* petsc = create_petsc_krylov_factory();
  test_2d_laplace_eqn(state, petsc, BICGSTAB_SOLVER, false);
}

void test_petsc_bicgstab_block_2d_laplace_eqn(void** state)
{
  krylov_factory_t* petsc = create_petsc_krylov_factory();
  test_2d_laplace_eqn(state, petsc, BICGSTAB_SOLVER, true);
}

void test_hypre_krylov_factory(void** state)
{
  krylov_factory_t* hypre = create_hypre_krylov_factory();
  test_krylov_factory(state, hypre);
}

void test_hypre_krylov_matrix(void** state)
{
  krylov_factory_t* hypre = create_hypre_krylov_factory();
  test_krylov_matrix(state, hypre);
}

void test_hypre_krylov_matrix_from_file(void** state)
{
  krylov_factory_t* hypre = create_hypre_krylov_factory();
  test_krylov_matrix_from_sherman1(state, hypre);
}

void test_hypre_krylov_vector(void** state)
{
  krylov_factory_t* hypre = create_hypre_krylov_factory();
  test_krylov_vector(state, hypre);
}

void test_hypre_krylov_vector_from_file(void** state)
{
  krylov_factory_t* hypre = create_hypre_krylov_factory();
  test_krylov_vector_from_sherman1_b(state, hypre);
}

void test_hypre_pcg_1d_laplace_eqn(void** state)
{
  krylov_factory_t* hypre = create_hypre_krylov_factory();
  test_1d_laplace_eqn(state, hypre, PCG_SOLVER);
}

void test_hypre_gmres_1d_laplace_eqn(void** state)
{
  krylov_factory_t* hypre = create_hypre_krylov_factory();
  test_1d_laplace_eqn(state, hypre, GMRES_SOLVER);
}

void test_hypre_bicgstab_1d_laplace_eqn(void** state)
{
  krylov_factory_t* hypre = create_hypre_krylov_factory();
  test_1d_laplace_eqn(state, hypre, BICGSTAB_SOLVER);
}

void test_hypre_sherman1(void** state)
{
  krylov_factory_t* hypre = create_hypre_krylov_factory();
  test_load_and_solve(state, hypre, 
                      CMAKE_CURRENT_SOURCE_DIR "/sherman1.mtx", 
                      CMAKE_CURRENT_SOURCE_DIR "/sherman1_b.mtx");
}

void test_hypre_pcg_2d_laplace_eqn(void** state)
{
  krylov_factory_t* hypre = create_hypre_krylov_factory();
  test_2d_laplace_eqn(state, hypre, PCG_SOLVER, false);
}

void test_hypre_pcg_block_2d_laplace_eqn(void** state)
{
  krylov_factory_t* hypre = create_hypre_krylov_factory();
  test_2d_laplace_eqn(state, hypre, PCG_SOLVER, true);
}

void test_hypre_gmres_2d_laplace_eqn(void** state)
{
  krylov_factory_t* hypre = create_hypre_krylov_factory();
  test_2d_laplace_eqn(state, hypre, GMRES_SOLVER, false);
}

void test_hypre_gmres_block_2d_laplace_eqn(void** state)
{
  krylov_factory_t* hypre = create_hypre_krylov_factory();
  test_2d_laplace_eqn(state, hypre, GMRES_SOLVER, true);
}

void test_hypre_bicgstab_2d_laplace_eqn(void** state)
{
  krylov_factory_t* hypre = create_hypre_krylov_factory();
  test_2d_laplace_eqn(state, hypre, BICGSTAB_SOLVER, false);
}

void test_hypre_bicgstab_block_2d_laplace_eqn(void** state)
{
  krylov_factory_t* hypre = create_hypre_krylov_factory();
  test_2d_laplace_eqn(state, hypre, BICGSTAB_SOLVER, true);
}
#endif

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_lis_krylov_factory),
    cmocka_unit_test(test_lis_krylov_matrix),
    cmocka_unit_test(test_lis_krylov_matrix_from_file),
    cmocka_unit_test(test_lis_krylov_vector),
    cmocka_unit_test(test_lis_krylov_vector_from_file),
    cmocka_unit_test(test_lis_pcg_1d_laplace_eqn),
    cmocka_unit_test(test_lis_gmres_1d_laplace_eqn),
    cmocka_unit_test(test_lis_bicgstab_1d_laplace_eqn),
    cmocka_unit_test(test_lis_sherman1),
    cmocka_unit_test(test_lis_pcg_2d_laplace_eqn),
    cmocka_unit_test(test_lis_pcg_block_2d_laplace_eqn),
    cmocka_unit_test(test_lis_gmres_2d_laplace_eqn),
    cmocka_unit_test(test_lis_gmres_block_2d_laplace_eqn),
    cmocka_unit_test(test_lis_bicgstab_2d_laplace_eqn),
    cmocka_unit_test(test_lis_bicgstab_block_2d_laplace_eqn),
#if POLYMEC_HAVE_SHARED_LIBS
    cmocka_unit_test(test_petsc_krylov_factory),
    cmocka_unit_test(test_petsc_krylov_matrix),
    cmocka_unit_test(test_petsc_krylov_matrix_from_file),
    cmocka_unit_test(test_petsc_krylov_vector),
    cmocka_unit_test(test_petsc_krylov_vector_from_file),
    cmocka_unit_test(test_petsc_pcg_1d_laplace_eqn),
    cmocka_unit_test(test_petsc_gmres_1d_laplace_eqn),
    cmocka_unit_test(test_petsc_bicgstab_1d_laplace_eqn),
    cmocka_unit_test(test_petsc_sherman1),
    cmocka_unit_test(test_petsc_pcg_2d_laplace_eqn),
    cmocka_unit_test(test_petsc_pcg_block_2d_laplace_eqn),
    cmocka_unit_test(test_petsc_gmres_2d_laplace_eqn),
    cmocka_unit_test(test_petsc_gmres_block_2d_laplace_eqn),
    cmocka_unit_test(test_petsc_bicgstab_2d_laplace_eqn),
    cmocka_unit_test(test_petsc_bicgstab_block_2d_laplace_eqn),
    cmocka_unit_test(test_hypre_krylov_factory),
    cmocka_unit_test(test_hypre_krylov_matrix),
    cmocka_unit_test(test_hypre_krylov_matrix_from_file),
    cmocka_unit_test(test_hypre_krylov_vector),
    cmocka_unit_test(test_hypre_krylov_vector_from_file),
    cmocka_unit_test(test_hypre_pcg_1d_laplace_eqn),
    cmocka_unit_test(test_hypre_gmres_1d_laplace_eqn),
    cmocka_unit_test(test_hypre_bicgstab_1d_laplace_eqn),
    cmocka_unit_test(test_hypre_sherman1),
    cmocka_unit_test(test_hypre_pcg_2d_laplace_eqn),
    cmocka_unit_test(test_hypre_pcg_block_2d_laplace_eqn),
    cmocka_unit_test(test_hypre_gmres_2d_laplace_eqn),
    cmocka_unit_test(test_hypre_gmres_block_2d_laplace_eqn),
    cmocka_unit_test(test_hypre_bicgstab_2d_laplace_eqn),
    cmocka_unit_test(test_hypre_bicgstab_block_2d_laplace_eqn)
#endif
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}

