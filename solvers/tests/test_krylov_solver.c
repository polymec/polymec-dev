// Copyright (c) 2012-2018, Jeffrey N. Johnson
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
#include "solvers/krylov_solver.h"

// This helper creates a sparsity pattern for a finite difference 
// discretization of the Laplacian operator. This pattern can be 
// used to create a 1D non-block representation, or a 2D block 
// representation of the laplacian operator.
static matrix_sparsity_t* create_laplacian_sparsity(MPI_Comm comm, int N_local)
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
    matrix_sparsity_t* sparsity = create_laplacian_sparsity(comm, 1000);

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

static void test_krylov_matrix_ops(void** state, krylov_factory_t* factory)
{
  if (factory != NULL)
  {
    // Create a serial matrix with 10 rows per process.
    int N = 10;
    matrix_sparsity_t* sparsity = create_laplacian_sparsity(MPI_COMM_SELF, 10);
    index_t* row_dist = matrix_sparsity_row_distribution(sparsity);
    index_t N_global = matrix_sparsity_num_global_rows(sparsity);
    krylov_matrix_t* A = krylov_factory_matrix(factory, sparsity);

    // Construct a test Laplacian to use with our operations.
    real_t h = 1.0 / N;
    int rpos = 0;
    index_t row;
    while (matrix_sparsity_next_row(sparsity, &rpos, &row))
    {
      real_t Aij[3];
      index_t rows[3], cols[3], num_cols;
      if (row == 0)
      {
        num_cols = 2; 
        rows[0] = 0; cols[0] = 0; Aij[0] = -2.0;
                     cols[1] = 1; Aij[1] = 1.0;
        krylov_matrix_set_values(A, 1, &num_cols, rows, cols, Aij);
      }
      else if (row == (N_global-1))
      {
        num_cols = 2;
        rows[0] = row; cols[0] = row-1; Aij[0] = 1.0;
                       cols[1] = row;   Aij[1] = -2.0;
        krylov_matrix_set_values(A, 1, &num_cols, rows, cols, Aij);
      }
      else
      {
        num_cols = 3;
        rows[0] = row; cols[0] = row-1; Aij[0] = 1.0;
        rows[1] = row; cols[1] = row;   Aij[1] = -2.0;
        rows[2] = row; cols[2] = row+1; Aij[2] = 1.0;
        krylov_matrix_set_values(A, 1, &num_cols, rows, cols, Aij);
      }
    }
    krylov_matrix_assemble(A);

    // Clone A and scale it.
    krylov_matrix_t* A1 = krylov_matrix_clone(A);
    krylov_matrix_scale(A1, 1.0/(h*h));
    rpos = 0;
    while (matrix_sparsity_next_row(sparsity, &rpos, &row))
    {
      real_t Aij[3];
      index_t rows[3], cols[3], num_cols;
      if (row == 0)
      {
        num_cols = 2; 
        rows[0] = 0; cols[0] = 0; cols[1] = 1;
        krylov_matrix_get_values(A1, 1, &num_cols, rows, cols, Aij);
        assert_true(reals_equal(Aij[0], -2.0/(h*h)));
        assert_true(reals_equal(Aij[1],  1.0/(h*h)));
      }
      else if (row == (N_global-1))
      {
        num_cols = 2;
        rows[0] = row; cols[0] = row-1; cols[1] = row;
        krylov_matrix_get_values(A1, 1, &num_cols, rows, cols, Aij);
        assert_true(reals_equal(Aij[0],  1.0/(h*h)));
        assert_true(reals_equal(Aij[1], -2.0/(h*h)));
      }
      else
      {
        num_cols = 3;
        rows[0] = row; cols[0] = row-1; cols[1] = row; cols[2] = row+1; 
        krylov_matrix_get_values(A1, 1, &num_cols, rows, cols, Aij);
        assert_true(reals_equal(Aij[0],  1.0/(h*h)));
        assert_true(reals_equal(Aij[1], -2.0/(h*h)));
        assert_true(reals_equal(Aij[2],  1.0/(h*h)));
      }
    }

    // Copy A1 <- A and diagonally scale.
    {
      krylov_matrix_copy(A, A1);
      krylov_vector_t* L = krylov_factory_vector(factory, MPI_COMM_SELF, row_dist);
      krylov_vector_t* R = krylov_vector_clone(L);
      real_t Li[10] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0}, 
      Ri[10] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};
      krylov_vector_copy_in(L, Li);
      krylov_vector_copy_in(R, Ri);
      krylov_matrix_diag_scale(A1, L, R);

      rpos = 0;
      while (matrix_sparsity_next_row(sparsity, &rpos, &row))
      {
        real_t Aij[3];
        index_t rows[3], cols[3], num_cols;
        if (row == 0)
        {
          num_cols = 2; 
          rows[0] = 0; cols[0] = 0; cols[1] = 1;
          krylov_matrix_get_values(A1, 1, &num_cols, rows, cols, Aij);
          assert_true(reals_equal(Aij[0], -Li[0]*2.0*Ri[0]));
          assert_true(reals_equal(Aij[1],  Li[0]*1.0*Ri[1]));
        }
        else if (row == (N_global-1))
        {
          num_cols = 2;
          rows[0] = row; cols[0] = row-1; cols[1] = row;
          krylov_matrix_get_values(A1, 1, &num_cols, rows, cols, Aij);
          assert_true(reals_equal(Aij[0], -Li[row]*2.0*Ri[cols[0]]));
          assert_true(reals_equal(Aij[1],  Li[row]*1.0*Ri[cols[1]]));
        }
        else
        {
          num_cols = 3;
          rows[0] = row; cols[0] = row-1; cols[1] = row; cols[2] = row+1; 
          krylov_matrix_get_values(A1, 1, &num_cols, rows, cols, Aij);
          assert_true(reals_equal(Aij[0],  Li[row]*1.0*Ri[cols[0]]));
          assert_true(reals_equal(Aij[1], -Li[row]*2.0*Ri[cols[1]]));
          assert_true(reals_equal(Aij[2],  Li[row]*1.0*Ri[cols[2]]));
        }
      }

      krylov_vector_free(L);
      krylov_vector_free(R);
    }

    // Copy A1 <- A and remove the diagonal.
    {
      krylov_matrix_copy(A, A1);
      krylov_matrix_add_identity(A1, 2.0);

      rpos = 0;
      while (matrix_sparsity_next_row(sparsity, &rpos, &row))
      {
        real_t Aij[3];
        index_t rows[3], cols[3], num_cols;
        if (row == 0)
        {
          num_cols = 2; 
          rows[0] = 0; cols[0] = 0; cols[1] = 1;
          krylov_matrix_get_values(A1, 1, &num_cols, rows, cols, Aij);
          assert_true(reals_equal(Aij[0], 0.0));
          assert_true(reals_equal(Aij[1], 1.0));
        }
        else if (row == (N_global-1))
        {
          num_cols = 2;
          rows[0] = row; cols[0] = row-1; cols[1] = row;
          krylov_matrix_get_values(A1, 1, &num_cols, rows, cols, Aij);
          assert_true(reals_equal(Aij[0], 1.0));
          assert_true(reals_equal(Aij[1], 0.0));
        }
        else
        {
          num_cols = 3;
          rows[0] = row; cols[0] = row-1; cols[1] = row; cols[2] = row+1;
          krylov_matrix_get_values(A1, 1, &num_cols, rows, cols, Aij);
          assert_true(reals_equal(Aij[0], 1.0));
          assert_true(reals_equal(Aij[1], 0.0));
          assert_true(reals_equal(Aij[2], 1.0));
        }
      }
    }

    // Clean up.
    krylov_matrix_free(A1);
    krylov_matrix_free(A);
    matrix_sparsity_free(sparsity);
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

static void test_krylov_vector_ops(void** state, krylov_factory_t* factory)
{
  if (factory != NULL)
  {
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
    matrix_sparsity_t* sparsity = create_laplacian_sparsity(comm, N);
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
      index_t rows[1], cols[3], num_cols;
      if (row == 0)
      {
        num_cols = 2; 
        rows[0] = 0; cols[0] = 0; Aij[0] = -2.0;
                     cols[1] = 1; Aij[1] = 1.0;
        krylov_matrix_set_values(A, 1, &num_cols, rows, cols, Aij);
        bi[0] = 0.0;
        krylov_vector_set_values(b, 1, &row, bi);
      }
      else if (row == (N_global-1))
      {
        num_cols = 2;
        rows[0] = row; cols[0] = row-1; Aij[0] = 1.0;
                       cols[1] = row;   Aij[1] = -2.0;
        krylov_matrix_set_values(A, 1, &num_cols, rows, cols, Aij);
        bi[0] = 1.0;
        krylov_vector_set_values(b, 1, &row, bi);
      }
      else
      {
        num_cols = 3;
        rows[0] = row; cols[0] = row-1; Aij[0] = 1.0;
                       cols[1] = row;   Aij[1] = -2.0;
                       cols[2] = row+1; Aij[2] = 1.0;
        krylov_matrix_set_values(A, 1, &num_cols, rows, cols, Aij);
        bi[0] = 0.0;
        krylov_vector_set_values(b, 1, &row, bi);
      }
    }
    krylov_matrix_assemble(A);
    krylov_vector_assemble(b);
    krylov_matrix_scale(A, 1.0/(h*h));
    krylov_vector_scale(b, -1.0/(h*h));

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
    matrix_sparsity_free(sparsity);
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

static void test_10x10_block(void** state, 
                             krylov_factory_t* factory)
{
  if (factory != NULL)
  {
    MPI_Comm comm = MPI_COMM_SELF;

    // Create a 10x10 matrix consisting of a single block.
    // Note that this matrix is actually the transpose of the block that 
    // will be inserted, since we assume column-major ordering matrices 
    // (for ease of interoperability with LAPACK).
    int N = 10;
    static real_t B[10*10] = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0,10.0,
                              11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0,20.0,
                              21.0,22.0,23.0,24.0,25.0,26.0,27.0,28.0,29.0,30.0,
                              31.0,32.0,33.0,34.0,35.0,36.0,37.0,38.0,39.0,40.0,
                              41.0,42.0,43.0,44.0,45.0,46.0,47.0,48.0,49.0,50.0,
                              51.0,52.0,53.0,54.0,55.0,56.0,57.0,58.0,59.0,60.0,
                              61.0,62.0,63.0,64.0,65.0,66.0,67.0,68.0,69.0,70.0,
                              71.0,72.0,73.0,74.0,75.0,76.0,77.0,78.0,79.0,80.0,
                              81.0,82.0,83.0,84.0,85.0,86.0,87.0,88.0,89.0,90.0,
                              91.0,92.0,93.0,94.0,95.0,96.0,97.0,98.0,99.0,100.0};
    index_t row_dist[2] = {0, 1};
    matrix_sparsity_t* sparsity = matrix_sparsity_new(comm, row_dist);
    matrix_sparsity_set_num_columns(sparsity, 0, 1);
    index_t* columns = matrix_sparsity_columns(sparsity, 0);
    columns[0] = 0;
    krylov_matrix_t* A = krylov_factory_block_matrix(factory, sparsity, N);
    matrix_sparsity_free(sparsity);

    // Shove the block into the matrix.
    krylov_matrix_set_block(A, 0, 0, B);
    krylov_matrix_assemble(A);
    krylov_matrix_free(A);
    krylov_factory_free(factory);
  }
}

static void test_2d_laplace_eqn(void** state, 
                                krylov_factory_t* factory,
                                iterative_solver_t solver_type)
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
    static real_t I[10*10] = {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                              0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                              0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                              0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                              0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
                              0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 
                              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 
                              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 
                              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 
                              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0};
    matrix_sparsity_t* sparsity = create_laplacian_sparsity(comm, N);
    index_t* block_row_dist = matrix_sparsity_row_distribution(sparsity);
    krylov_matrix_t* A = krylov_factory_block_matrix(factory, sparsity, N);
    index_t N_global = krylov_matrix_num_global_rows(A);

    // Create a RHS vector.
    index_t row_dist[nprocs+1];
    for (int p = 0; p <= nprocs; ++p)
      row_dist[p] = N * block_row_dist[p];
    krylov_vector_t* b = krylov_factory_vector(factory, comm, row_dist);

    // Construct the 2D Laplacian operator and RHS.
    real_t h = 1.0 / N;
    {
      int rpos = 0;
      index_t row;
      while (matrix_sparsity_next_row(sparsity, &rpos, &row))
      {
        // Fill the row of the matrix.
        if (row == 0)
        {
          krylov_matrix_set_block(A, 0, 0, T);
          krylov_matrix_set_block(A, 0, 1, I);
        }
        else if (row == (nprocs*N-1))
        {
          krylov_matrix_set_block(A, nprocs*N-1, nprocs*N-2, I);
          krylov_matrix_set_block(A, nprocs*N-1, nprocs*N-1, T);
        }
        else
        {
          krylov_matrix_set_block(A, row, row-1, I);
          krylov_matrix_set_block(A, row, row, T);
          krylov_matrix_set_block(A, row, row+1, I);
        }
      }
    }
    krylov_matrix_assemble(A);

    // Fill the rows of the RHS vector.
    for (index_t row = row_dist[rank]; row < row_dist[rank+1]; ++row)
    {
      real_t bi = 0.0;
      if (row >= N_global - N)
        bi = 1.0;
      else if (((row+1) % 10) == 0)
        bi = 1.0;
      krylov_vector_set_values(b, 1, &row, &bi);
    }
    krylov_vector_assemble(b);

    krylov_matrix_scale(A, 1.0/(h*h));
    krylov_vector_scale(b, -1.0/(h*h));

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
    matrix_sparsity_free(sparsity);
  }
}

static void test_petsc_krylov_factory(void** state)
{
  krylov_factory_t* petsc = create_petsc_krylov_factory();
  test_krylov_factory(state, petsc);
}

static void test_petsc_krylov_matrix(void** state)
{
  krylov_factory_t* petsc = create_petsc_krylov_factory();
  test_krylov_matrix(state, petsc);
}

static void test_petsc_krylov_matrix_from_file(void** state)
{
  krylov_factory_t* petsc = create_petsc_krylov_factory();
  test_krylov_matrix_from_sherman1(state, petsc);
}

static void test_petsc_krylov_matrix_ops(void** state)
{
  krylov_factory_t* petsc = create_petsc_krylov_factory();
  test_krylov_matrix_ops(state, petsc);
}

static void test_petsc_krylov_vector(void** state)
{
  krylov_factory_t* petsc = create_petsc_krylov_factory();
  test_krylov_vector(state, petsc);
}

static void test_petsc_krylov_vector_from_file(void** state)
{
  krylov_factory_t* petsc = create_petsc_krylov_factory();
  test_krylov_vector_from_sherman1_b(state, petsc);
}

static void test_petsc_krylov_vector_ops(void** state)
{
  krylov_factory_t* petsc = create_petsc_krylov_factory();
  test_krylov_vector_ops(state, petsc);
}

static void test_petsc_pcg_1d_laplace_eqn(void** state)
{
  krylov_factory_t* petsc = create_petsc_krylov_factory();
  test_1d_laplace_eqn(state, petsc, PCG_SOLVER);
}

static void test_petsc_gmres_1d_laplace_eqn(void** state)
{
  krylov_factory_t* petsc = create_petsc_krylov_factory();
  test_1d_laplace_eqn(state, petsc, GMRES_SOLVER);
}

static void test_petsc_bicgstab_1d_laplace_eqn(void** state)
{
  krylov_factory_t* petsc = create_petsc_krylov_factory();
  test_1d_laplace_eqn(state, petsc, BICGSTAB_SOLVER);
}

static void test_petsc_sherman1(void** state)
{
  krylov_factory_t* petsc = create_petsc_krylov_factory();
  test_load_and_solve(state, petsc, 
                      CMAKE_CURRENT_SOURCE_DIR "/sherman1.mtx", 
                      CMAKE_CURRENT_SOURCE_DIR "/sherman1_b.mtx");
}

static void test_petsc_10x10_block(void** state)
{
  krylov_factory_t* petsc = create_petsc_krylov_factory();
  test_10x10_block(state, petsc);
}

static void test_petsc_pcg_2d_laplace_eqn(void** state)
{
  krylov_factory_t* petsc = create_petsc_krylov_factory();
  test_2d_laplace_eqn(state, petsc, PCG_SOLVER);
}

static void test_petsc_gmres_2d_laplace_eqn(void** state)
{
  krylov_factory_t* petsc = create_petsc_krylov_factory();
  test_2d_laplace_eqn(state, petsc, GMRES_SOLVER);
}

static void test_petsc_bicgstab_2d_laplace_eqn(void** state)
{
  krylov_factory_t* petsc = create_petsc_krylov_factory();
  test_2d_laplace_eqn(state, petsc, BICGSTAB_SOLVER);
}

static void test_hypre_krylov_factory(void** state)
{
  krylov_factory_t* hypre = create_hypre_krylov_factory();
  test_krylov_factory(state, hypre);
}

static void test_hypre_krylov_matrix(void** state)
{
  krylov_factory_t* hypre = create_hypre_krylov_factory();
  test_krylov_matrix(state, hypre);
}

static void test_hypre_krylov_matrix_from_file(void** state)
{
  krylov_factory_t* hypre = create_hypre_krylov_factory();
  test_krylov_matrix_from_sherman1(state, hypre);
}

static void test_hypre_krylov_matrix_ops(void** state)
{
  krylov_factory_t* hypre = create_hypre_krylov_factory();
  test_krylov_matrix_ops(state, hypre);
}

static void test_hypre_krylov_vector(void** state)
{
  krylov_factory_t* hypre = create_hypre_krylov_factory();
  test_krylov_vector(state, hypre);
}

static void test_hypre_krylov_vector_from_file(void** state)
{
  krylov_factory_t* hypre = create_hypre_krylov_factory();
  test_krylov_vector_from_sherman1_b(state, hypre);
}

static void test_hypre_krylov_vector_ops(void** state)
{
  krylov_factory_t* hypre = create_hypre_krylov_factory();
  test_krylov_vector_ops(state, hypre);
}

static void test_hypre_pcg_1d_laplace_eqn(void** state)
{
  krylov_factory_t* hypre = create_hypre_krylov_factory();
  test_1d_laplace_eqn(state, hypre, PCG_SOLVER);
}

static void test_hypre_gmres_1d_laplace_eqn(void** state)
{
  krylov_factory_t* hypre = create_hypre_krylov_factory();
  test_1d_laplace_eqn(state, hypre, GMRES_SOLVER);
}

static void test_hypre_bicgstab_1d_laplace_eqn(void** state)
{
  krylov_factory_t* hypre = create_hypre_krylov_factory();
  test_1d_laplace_eqn(state, hypre, BICGSTAB_SOLVER);
}

static void test_hypre_sherman1(void** state)
{
  krylov_factory_t* hypre = create_hypre_krylov_factory();
  test_load_and_solve(state, hypre, 
                      CMAKE_CURRENT_SOURCE_DIR "/sherman1.mtx", 
                      CMAKE_CURRENT_SOURCE_DIR "/sherman1_b.mtx");
}

static void test_hypre_10x10_block(void** state)
{
  krylov_factory_t* hypre = create_hypre_krylov_factory();
  test_10x10_block(state, hypre);
}


static void test_hypre_pcg_2d_laplace_eqn(void** state)
{
  krylov_factory_t* hypre = create_hypre_krylov_factory();
  test_2d_laplace_eqn(state, hypre, PCG_SOLVER);
}

static void test_hypre_gmres_2d_laplace_eqn(void** state)
{
  krylov_factory_t* hypre = create_hypre_krylov_factory();
  test_2d_laplace_eqn(state, hypre, GMRES_SOLVER);
}

static void test_hypre_bicgstab_2d_laplace_eqn(void** state)
{
  krylov_factory_t* hypre = create_hypre_krylov_factory();
  test_2d_laplace_eqn(state, hypre, BICGSTAB_SOLVER);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_petsc_krylov_factory),
    cmocka_unit_test(test_petsc_krylov_matrix),
    cmocka_unit_test(test_petsc_krylov_matrix_from_file),
    cmocka_unit_test(test_petsc_krylov_matrix_ops),
    cmocka_unit_test(test_petsc_krylov_vector),
    cmocka_unit_test(test_petsc_krylov_vector_from_file),
    cmocka_unit_test(test_petsc_krylov_vector_ops),
    cmocka_unit_test(test_petsc_pcg_1d_laplace_eqn),
    cmocka_unit_test(test_petsc_gmres_1d_laplace_eqn),
    cmocka_unit_test(test_petsc_bicgstab_1d_laplace_eqn),
    cmocka_unit_test(test_petsc_sherman1),
    cmocka_unit_test(test_petsc_10x10_block),
    cmocka_unit_test(test_petsc_pcg_2d_laplace_eqn),
    cmocka_unit_test(test_petsc_gmres_2d_laplace_eqn),
    cmocka_unit_test(test_petsc_bicgstab_2d_laplace_eqn),
    cmocka_unit_test(test_hypre_krylov_factory),
    cmocka_unit_test(test_hypre_krylov_matrix),
    cmocka_unit_test(test_hypre_krylov_matrix_from_file),
    cmocka_unit_test(test_hypre_krylov_matrix_ops),
    cmocka_unit_test(test_hypre_krylov_vector),
    cmocka_unit_test(test_hypre_krylov_vector_from_file),
    cmocka_unit_test(test_hypre_krylov_vector_ops),
    cmocka_unit_test(test_hypre_pcg_1d_laplace_eqn),
    cmocka_unit_test(test_hypre_gmres_1d_laplace_eqn),
    cmocka_unit_test(test_hypre_bicgstab_1d_laplace_eqn),
    cmocka_unit_test(test_hypre_sherman1),
    cmocka_unit_test(test_hypre_10x10_block),
    cmocka_unit_test(test_hypre_pcg_2d_laplace_eqn),
    cmocka_unit_test(test_hypre_gmres_2d_laplace_eqn),
    cmocka_unit_test(test_hypre_bicgstab_2d_laplace_eqn)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}

