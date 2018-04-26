// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_KRYLOV_SOLVER_H
#define POLYMEC_KRYLOV_SOLVER_H

#include "core/polymec.h"
#include "core/unordered_map.h"
#include "solvers/matrix_sparsity.h"

// The Krylov solver interface is an abstract interface that can be used with 
// third-party parallel sparse linear solvers. The linkage to these solvers is 
// achieved with dynamic loading, so this interface is only useful on platforms
// that support dynamic loading.

// Objects of this type store distributed vectors for use with Krylov solvers.
typedef struct krylov_vector_t krylov_vector_t;

// Objects of this type store distributed matrices for use with Krylov solvers.
typedef struct krylov_matrix_t krylov_matrix_t;

// Objects of this type solve distributed linear systems A*x = b using 
// preconditioned Krylov subspace methods.
typedef struct krylov_solver_t krylov_solver_t;

// Objects of this type are preconditioners for the krylov_solvers.
typedef struct krylov_pc_t krylov_pc_t;

// Objects of this type construct Krylov matrices, vectors, solvers, and preconditioners.
// A Krylov factory exposes capabilities within a given library (PETSc, HYPRE, etc). The 
// factory must continue to exist as long as any of the solvers, matrices, and vectors it 
// has created.
typedef struct krylov_factory_t krylov_factory_t;

//------------------------------------------------------------------------
//        Machinery for implementing third-party Krylov solvers.
//------------------------------------------------------------------------

// This virtual table must be filled out for any subclass of krylov_factory.
typedef struct 
{
  krylov_solver_t* (*pcg_solver)(void* context, MPI_Comm comm);
  krylov_solver_t* (*gmres_solver)(void* context, MPI_Comm comm, int krylov_dimension);
  krylov_solver_t* (*bicgstab_solver)(void* context, MPI_Comm comm);
  krylov_solver_t* (*special_solver)(void* context, MPI_Comm comm, const char* solver_name, string_string_unordered_map_t* options);
  krylov_pc_t* (*preconditioner)(void* context, MPI_Comm comm, const char* pc_name, string_string_unordered_map_t* options);
  krylov_matrix_t* (*matrix)(void* context, matrix_sparsity_t* sparsity);
  krylov_matrix_t* (*block_matrix)(void* context, matrix_sparsity_t* sparsity, size_t block_size);
  krylov_matrix_t* (*var_block_matrix)(void* context, matrix_sparsity_t* sparsity, size_t* block_sizes);
  krylov_vector_t* (*vector)(void* context, MPI_Comm comm, index_t* row_dist);
  void (*dtor)(void* context);
} krylov_factory_vtable;

// This constructor should be called with a context pointer and a virtual 
// table to create an instance of a krylov_factory subclass.
krylov_factory_t* krylov_factory_new(const char* name,
                                     void* context,
                                     krylov_factory_vtable vtable);

// This virtual table must be filled out for any subclass of krylov_solver.
typedef struct
{
  void (*set_tolerances)(void* context, real_t rel_tol, real_t abs_tol, real_t div_tol);
  void (*set_max_iterations)(void* context, int max_iters);
  void (*set_operator)(void* context, void* op);
  void (*set_preconditioner)(void* context, void* pc);
  bool (*solve)(void* context, void* b, void* x, real_t* res_norm, int* num_iters);
  bool (*solve_scaled)(void* context, void* b, void* s1, void* s2, void* x, real_t* res_norm, int* num_iters);
  void (*dtor)(void* context);
} krylov_solver_vtable;

// This constructor should be called with a context pointer and a virtual 
// table to create an instance of a krylov_solver subclass.
krylov_solver_t* krylov_solver_new(const char* name,
                                   void* context,
                                   krylov_solver_vtable vtable);

// This virtual table must be filled out for any subclass of krylov_pc.
typedef struct
{
  void (*dtor)(void* context);
} krylov_pc_vtable;

// This constructor should be called with a context pointer and a virtual table 
// to create an instance of a krylov_pc subclass.
krylov_pc_t* krylov_pc_new(const char* name,
                           void* context,
                           krylov_pc_vtable vtable);

// This virtual table must be filled out for any subclass of krylov_matrix.
typedef struct
{
  size_t (*block_size)(void* context, index_t block_row);
  void* (*clone)(void* context);
  void (*copy)(void* context, void* copy);
  void (*redistribute)(void* context, MPI_Comm);
  void (*zero)(void* context);
  void (*scale)(void* context, real_t scale_factor);
  void (*diag_scale)(void* context, void* L, void* R);
  void (*add_identity)(void* context, real_t scale_factor);
  void (*add_diagonal)(void* context, void* D);
  void (*set_diagonal)(void* context, void* D);
  void (*matvec)(void* context, void* x, bool transpose, void* y);
  void (*set_values)(void* context, size_t num_rows, size_t* num_columns, index_t* rows, index_t* columns, real_t* values);
  void (*add_values)(void* context, size_t num_rows, size_t* num_columns, index_t* rows, index_t* columns, real_t* values);
  void (*get_values)(void* context, size_t num_rows, size_t* num_columns, index_t* rows, index_t* columns, real_t* values);
  void (*set_blocks)(void* context, size_t num_blocks, index_t* block_rows, index_t* block_columns, real_t* block_values);
  void (*add_blocks)(void* context, size_t num_blocks, index_t* block_rows, index_t* block_columns, real_t* block_values);
  void (*get_blocks)(void* context, size_t num_blocks, index_t* block_rows, index_t* block_columns, real_t* block_values);
  void (*assemble)(void* context);
  void (*fprintf)(void* context, FILE* stream);
  void (*dtor)(void* context);
} krylov_matrix_vtable;

// This constructor should be called with a context pointer and a virtual 
// table to create an instance of a krylov_matrix subclass.
krylov_matrix_t* krylov_matrix_new(void* context,
                                   krylov_matrix_vtable vtable,
                                   MPI_Comm comm,
                                   size_t num_local_rows,
                                   size_t num_global_rows);

// This virtual table must be filled out for any subclass of krylov_matrix.
typedef struct
{
  void* (*clone)(void* context);
  void (*copy)(void* context, void* copy);
  void (*zero)(void* context);
  void (*set_value)(void* context, real_t value);
  void (*scale)(void* context, real_t scale_factor);
  void (*diag_scale)(void* context, void* D);
  void (*set_values)(void* context, size_t num_values, index_t* indices, real_t* values);
  void (*add_values)(void* context, size_t num_values, index_t* indices, real_t* values);
  void (*get_values)(void* context, size_t num_values, index_t* indices, real_t* values);
  void (*copy_in)(void* context, real_t* local_values);
  void (*copy_out)(void* context, real_t* local_values);
  real_t (*dot)(void* context, void* w);
  real_t (*norm)(void* context, int p);
  real_t (*w2_norm)(void* context, void* w);
  real_t (*wrms_norm)(void* context, void* w);
  void (*assemble)(void* context);
  void (*fprintf)(void* context, FILE* stream);
  void (*dtor)(void* context);
} krylov_vector_vtable;

// This constructor should be called with a context pointer and a virtual 
// table to create an instance of a krylov_vector subclass.
krylov_vector_t* krylov_vector_new(void* context,
                                   krylov_vector_vtable vtable,
                                   size_t local_size,
                                   size_t global_size);

//------------------------------------------------------------------------
//                  Bundled Krylov factories 
//------------------------------------------------------------------------

// This creates a PETSc-based Krylov factory that can be used for constructing
// matrices, vectors, solvers. The factory attempts to load the dynamic library
// located at petsc_library. If no such library can be found or loaded, 
// this function returns NULL. Use the use_64_bit_indices flag to indicate 
// whether this PETSc library was built using --with-64-bit-indices (since 
// we can't infer this using dynamic loading).
krylov_factory_t* petsc_krylov_factory(const char* petsc_library,
                                       bool use_64_bit_indices);

// This creates a HYPRE-based Krylov factory that can be used for constructing
// matrices, vectors, solvers, using the dynamic library located at 
// hypre_library. Use the use_64_bit_indices flag to indicate 
// whether this PETSc library was built using --with-64-bit-indices (since 
// we can't infer this using dynamic loading).
krylov_factory_t* hypre_krylov_factory(const char* hypre_library,
                                       bool use_64_bit_indices);

//------------------------------------------------------------------------
//                    Krylov factory interface
//------------------------------------------------------------------------

// Destroys an existing krylov_factory.
void krylov_factory_free(krylov_factory_t* factory);

// Returns an internal string containing the name of the given Krylov 
// factory implementation.
char* krylov_factory_name(krylov_factory_t* factory);

// Constructs a (square) Krylov sparse matrix with the sparsity pattern given 
// by an adjacency graph.
krylov_matrix_t* krylov_factory_matrix(krylov_factory_t* factory, 
                                       matrix_sparsity_t* sparsity);

// Constructs a (square) Krylov sparse block matrix with the sparsity pattern 
// given by an adjacency graph, and the given block size.
krylov_matrix_t* krylov_factory_block_matrix(krylov_factory_t* factory, 
                                             matrix_sparsity_t* sparsity,
                                             size_t block_size);

// Constructs a (square) Krylov sparse block matrix with the sparsity pattern 
// given by an adjacency graph, and different block sizes for each row.
krylov_matrix_t* krylov_factory_var_block_matrix(krylov_factory_t* factory, 
                                                 matrix_sparsity_t* sparsity,
                                                 size_t* block_sizes);

// Reads a matrix into memory from the given file (assuming it is the 
// in a supported file format), distributing it over the processes
// in the given communicator using a naive partitioning. Only the Matrix Market 
// file format is currently supported.
krylov_matrix_t* krylov_factory_matrix_from_file(krylov_factory_t* factory,
                                                 MPI_Comm comm,
                                                 const char* filename);

// Constructs a vector on the given communicator with its local and global 
// dimensions defined by the given row distribution array (which is 
// nprocs + 1 in length). row_distribution[p] gives the global index of the 
// first row stored on process p, and the number of rows local to that 
// process is given by row_distribution[p+1] - row_distribution[p].
krylov_vector_t* krylov_factory_vector(krylov_factory_t* factory,
                                       MPI_Comm comm,
                                       index_t* row_distribution);

// Reads a vector into memory from the given file (assuming it is the 
// in a supported file format), distributing it over the processes
// in the given communicator using a naive partitioning. Only the Matrix Market 
// file format is currently supported.
krylov_vector_t* krylov_factory_vector_from_file(krylov_factory_t* factory,
                                                 MPI_Comm comm,
                                                 const char* filename);

// Constructs a preconditioned conjugate gradient (PCG) Krylov solver. Keep in 
// mind that this method can only be used for systems having symmetric, 
// positive-definite matrices.
krylov_solver_t* krylov_factory_pcg_solver(krylov_factory_t* factory,
                                           MPI_Comm comm);
                                             
// Constructs a GMRES Krylov solver with the given Krylov subspace dimension.
krylov_solver_t* krylov_factory_gmres_solver(krylov_factory_t* factory,
                                             MPI_Comm comm,
                                             int krylov_dimension);
                                             
// Constructs a stabilized bi-conjugate gradient Krylov solver.
krylov_solver_t* krylov_factory_bicgstab_solver(krylov_factory_t* factory,
                                                MPI_Comm comm);

// Constructs a special Krylov solver supported by the backend in use, 
// identified by its name, with options specified in string key-value 
// pairs. If the solver with the given name is not available, this 
// function returns NULL.
krylov_solver_t* krylov_factory_special_solver(krylov_factory_t* factory,
                                               MPI_Comm comm,
                                               const char* solver_name,
                                               string_string_unordered_map_t* options);

// Constructs a preconditioner supported by the backend in use, 
// identified by its name, with options specified in string key-value
// pairs. If the preconditioner with the given name is not available, this 
// function returns NULL.
krylov_pc_t* krylov_factory_preconditioner(krylov_factory_t* factory,
                                           MPI_Comm comm,
                                           const char* pc_name,
                                           string_string_unordered_map_t* options);

//------------------------------------------------------------------------
//                      Krylov solver interface
//------------------------------------------------------------------------

// Frees a solver.
void krylov_solver_free(krylov_solver_t* solver);

// Returns an internal string containing the name of the solver.
char* krylov_solver_name(krylov_solver_t* solver);

// Returns a pointer to the underlying solver implementation. You can use 
// this if you have explicitly linked your program to the library providing
// the implementation.
void* krylov_solver_impl(krylov_solver_t* solver);

// Sets the tolerances for the residual norm (norm_tolerance) = |R|, where 
// R = A*x - b.
void krylov_solver_set_tolerances(krylov_solver_t* solver, 
                                  real_t relative_tolerance,
                                  real_t absolute_tolerance,
                                  real_t divergence_tolerance);

// Sets the maximum number of iterations for a linear solve.
void krylov_solver_set_max_iterations(krylov_solver_t* solver, 
                                      int max_iterations);

// Sets the operator A in the linear system to be solved.
void krylov_solver_set_operator(krylov_solver_t* solver, 
                                krylov_matrix_t* op);

// Sets the preconditioner for the linear solver. The solver assumes 
// responsibility for destroying the preconditioner.
void krylov_solver_set_preconditioner(krylov_solver_t* solver,
                                      krylov_pc_t* preconditioner);

// Returns an internal pointer to the preconditioner for the linear solver, 
// or NULL if no preconditioner is set.
krylov_pc_t* krylov_solver_preconditioner(krylov_solver_t* solver);

// Solves the linear system of equations A * x = b, storing the solution
// in the vector x. Returns true if the solution was obtained, false if not.
// The operator matrix A must be set using krylov_solver_set_operator before
// this function is called. The residual norm ||A * x - b|| will be stored in 
// *residual_norm, and the number of linear iterations will be stored in 
// *num_iterations upon success. If the solve is unsuccessful, the vector x 
// will NOT contain the solution unless the residual was reduced.
bool krylov_solver_solve(krylov_solver_t* solver, 
                         krylov_vector_t* b, 
                         krylov_vector_t* x, 
                         real_t* residual_norm, 
                         int* num_iterations);

// Solves the scaled linear system of equations 
// (s1 * A * s2_inv) * (s2 * x) = (s1 * b), storing the unscaled solution x
// in the vector x. Returns true if the solution was obtained, false if not.
// s1 and s2 are diagonal matrices represented by vectors. Either or both of 
// them may be set to NULL, in which case they are assumed to be the unscaled 
// identity matrix. The operator matrix A must be set using 
// krylov_solver_set_operator before this function is called. The residual norm 
// ||s1 * (b - A * x)||_2 will be stored in *residual_norm, and the number of 
// linear iterations will be stored in *num_iterations upon success. If the 
// solve is unsuccessful, the vector x will NOT contain the solution unless the 
// residual was reduced. 
bool krylov_solver_solve_scaled(krylov_solver_t* solver, 
                                krylov_vector_t* b, 
                                krylov_vector_t* s1,
                                krylov_vector_t* s2,
                                krylov_vector_t* x, 
                                real_t* residual_norm, 
                                int* num_iterations);

//------------------------------------------------------------------------
//                      Krylov preconditioner interface
//------------------------------------------------------------------------

// Frees a preconditioner.
void krylov_pc_free(krylov_pc_t* preconditioner);

// Returns an internal string containing the name of the preconditioner.
char* krylov_pc_name(krylov_solver_t* preconditioner);

//------------------------------------------------------------------------
//                      Krylov matrix interface
//------------------------------------------------------------------------

// Frees a matrix.
void krylov_matrix_free(krylov_matrix_t* A);

// Returns the communicator on which the matrix is defined.
MPI_Comm krylov_matrix_comm(krylov_matrix_t* A);

// Returns the size of the given block row. If the matrix doesn't use block 
// rows, the block size is 1.
size_t krylov_matrix_block_size(krylov_matrix_t* A, index_t block_row);

// Creates and returns a deep copy of a matrix.
krylov_matrix_t* krylov_matrix_clone(krylov_matrix_t* A);

// Copies the contents of the matrix A to those of copy.
void krylov_matrix_copy(krylov_matrix_t* A, krylov_matrix_t* copy);

// Returns a newly created matrix with the same data, redistributed from the 
// existing matrix's communicator to the given one.
krylov_matrix_t* krylov_matrix_redistribute(krylov_matrix_t* A, 
                                            MPI_Comm comm);

// Returns a pointer to the underlying matrix implementation. You can use 
// this if you have explicitly linked your program to the library providing
// the implementation.
void* krylov_matrix_impl(krylov_matrix_t* A);

// Returns the number of locally stored rows in the matrix.
size_t krylov_matrix_num_local_rows(krylov_matrix_t* A);

// Returns the number of globally stored rows in the matrix.
size_t krylov_matrix_num_global_rows(krylov_matrix_t* A);

// Zeros all of the entries in the given matrix.
// This is collective and must be called by all processes.
void krylov_matrix_zero(krylov_matrix_t* A);

// Multiplies this matrix by a scale factor.
// This is collective and must be called by all processes.
void krylov_matrix_scale(krylov_matrix_t* A,
                         real_t scale_factor);

// Left-multiplies this matrix by a left diagonal matrix L, and 
// right-multiples by a diagonal matrix R. L and R are represented 
// by vectors. If L or R is NULL, it is assumed to be the identity 
// matrix. This is collective and must be called by all processes.
void krylov_matrix_diag_scale(krylov_matrix_t* A,
                              krylov_vector_t* L, 
                              krylov_vector_t* R);

// Adds a scaled identity matrix to this one.
// This is collective and must be called by all processes.
void krylov_matrix_add_identity(krylov_matrix_t* A,
                                real_t scale_factor);

// Adds the entries of the given vector (D) to the diagonal of the given matrix (A).
// This is collective and must be called by all processes.
void krylov_matrix_add_diagonal(krylov_matrix_t* A,
                                krylov_vector_t* D);

// Sets the diagonal entries in the given matrix (A) to those in the given vector (D).
// This is collective and must be called by all processes.
void krylov_matrix_set_diagonal(krylov_matrix_t* A,
                                krylov_vector_t* D);

// Computes the matrix-vector product A*x, storing the result in y. If transpose
// is set to true, then x is left-multiplied by the transpose of A--otherwise
// it is left-multiplied by A.
void krylov_matrix_matvec(krylov_matrix_t* A,
                          krylov_vector_t* x,
                          bool transpose,
                          krylov_vector_t* y);

// Sets the values of the elements in the matrix identified by the given 
// (globally-indexed) rows and columns. The rows and columns being set must 
// exist on the local process. num_columns and rows are arrays of size num_rows, 
// and contain the numbers of columns and the row indices for each row. columns
// is an array containing the indices of all the columns, in row-major, column-minor
// order. values is an array with values corresponding to the entries in the columns 
// array.
void krylov_matrix_set_values(krylov_matrix_t* A,
                              size_t num_rows,
                              size_t* num_columns,
                              index_t* rows, 
                              index_t* columns,
                              real_t* values);
                              
// Adds the given values of the elements to those in the matrix.
// The values are identified by the given (globally-indexed) rows and columns.
// These rows and columns must exist on the local process. num_columns and rows are 
// arrays of size num_rows, and contain the numbers of columns and the row indices for 
// each row. columns is an array containing the indices of all the columns, in row-major, 
// column-minor order. values is an array with values corresponding to the entries in the 
// columns array.
void krylov_matrix_add_values(krylov_matrix_t* A,
                              size_t num_rows,
                              size_t* num_columns,
                              index_t* rows, 
                              index_t* columns,
                              real_t* values);
                              
// Retrieves the values of the elements in the matrix identified by the 
// given (globally-indexed) rows and columns, storing them in the values array.
void krylov_matrix_get_values(krylov_matrix_t* A,
                              size_t num_rows,
                              size_t* num_columns,
                              index_t* rows, 
                              index_t* columns,
                              real_t* values);

// Sets the values of the blocks in the matrix identified by the given 
// (globally-indexed) block rows and columns. The block rows and columns being 
// set must exist on the local process. block_rows and block_columns are arrays 
// of size num_blocks, and contain the indices of the blocks being set. 
// block_values is an array of size block_size*block_size*num_blocks whose data
// consists of an array of num_blocks column-major-ordered blocks of data.
void krylov_matrix_set_blocks(krylov_matrix_t* A,
                              size_t num_blocks,
                              index_t* block_rows, 
                              index_t* block_columns,
                              real_t* block_values);
                              
// Adds in the values of the blocks in the matrix identified by the given 
// (globally-indexed) block rows and columns. The block rows and columns being 
// set must exist on the local process. block_rows and block_columns are arrays 
// of size num_blocks, and contain the indices of the blocks being set. 
// block_values is an array of size block_size*block_size*num_blocks whose data
// consists of an array of num_blocks column-major-ordered blocks of data.
void krylov_matrix_add_blocks(krylov_matrix_t* A,
                              size_t num_blocks,
                              index_t* block_rows, 
                              index_t* block_columns,
                              real_t* block_values);
                              
// Retrieves the blocks of values in the matrix identified by the given 
// (globally-indexed) block rows and columns, storing them in the block_values 
// array (as an array of num_blocks sets of block_size*block_size 
// column-major-ordered values). These blocks must be stored on the local 
// process.  
void krylov_matrix_get_blocks(krylov_matrix_t* A,
                              size_t num_blocks,
                              index_t* block_rows, 
                              index_t* block_columns,
                              real_t* block_values);

// Set the values of the block in the matrix identified by the given 
// (globally-indexed) block row and column, which must exist on the local 
// process. The block_values array contains block_size*block_size 
// column-major-ordered elements.
void krylov_matrix_set_block(krylov_matrix_t* A,
                             index_t block_row, 
                             index_t block_column,
                             real_t* block_values);
                              
// Adds in the values of the block in the matrix identified by the given 
// (globally-indexed) block row and column, which must exist on the local 
// process. The block_values array contains block_size*block_size 
// column-major-ordered elements.
void krylov_matrix_add_block(krylov_matrix_t* A,
                             index_t block_row, 
                             index_t block_column,
                             real_t* block_values);
                              
// Retrieves the block of values in the matrix identified by the given 
// (globally-indexed) block row and column, storing it in the block_values 
// array (as an array of block_size*block_size column-major-ordered values). 
// The block must be stored on the local process.  
void krylov_matrix_get_block(krylov_matrix_t* A,
                             index_t block_row, 
                             index_t block_column,
                             real_t* block_values);

// Assembles added/inserted values into the matrix, allowing all the processes
// a consistent representation of the matrix. This should be called after calls
// to krylov_matrix_set_values/krylov_matrix_add_values and 
// krylov_matrix_set_blocks/krylov_matrix_add_blocks, and should be 
// placed in between sets and adds.
void krylov_matrix_assemble(krylov_matrix_t* A);

// Writes a text representation of the matrix (or portion stored on the local
// MPI process) to the given stream.
void krylov_matrix_fprintf(krylov_matrix_t* A,
                           FILE* stream);

//------------------------------------------------------------------------
//                      Krylov vector interface
//------------------------------------------------------------------------

// Frees a vector.
void krylov_vector_free(krylov_vector_t* v);

// Creates and returns a deep copy of a vector.
krylov_vector_t* krylov_vector_clone(krylov_vector_t* v);

// Copies the contents of the vector v to those of copy.
void krylov_vector_copy(krylov_vector_t* v, krylov_vector_t* copy);

// Returns a pointer to the underlying vector implementation. You can use 
// this if you have explicitly linked your program to the library providing
// the implementation.
void* krylov_vector_impl(krylov_vector_t* v);

// Returns the locally-stored size (dimension) of the vector.
size_t krylov_vector_local_size(krylov_vector_t* v);

// Returns the global size (dimension) of the vector.
size_t krylov_vector_global_size(krylov_vector_t* v);

// Zeros all of the entries in the given vector.
// This is collective, and must be called by all MPI processes.
void krylov_vector_zero(krylov_vector_t* v);

// Sets all of the entries in the given vector to the given value.
// This is collective, and must be called by all MPI processes.
void krylov_vector_set_value(krylov_vector_t* v,
                             real_t value);

// Scales the vector by the given factor.
// This is collective, and must be called by all MPI processes.
void krylov_vector_scale(krylov_vector_t* v,
                         real_t scale_factor);

// Multiplies the vector by the diagonal matrix D, represented by a vector.
// This is collective, and must be called by all MPI processes.
void krylov_vector_diag_scale(krylov_vector_t* v,
                              krylov_vector_t* D);

// Sets the values of the elements in the vector identified by the given 
// (globally-indexed) indices.
void krylov_vector_set_values(krylov_vector_t* v,
                              size_t num_values,
                              index_t* indices,
                              real_t* values);
                              
// Adds the given values of the elements to those in the vector.
// The values are identified by the given (globally-indexed) indices.
void krylov_vector_add_values(krylov_vector_t* v,
                              size_t num_values,
                              index_t* indices,
                              real_t* values);

// Copies the data in the given array of local values into this vector. The 
// array is assumed to contain data for all rows of the vector stored locally.
void krylov_vector_copy_in(krylov_vector_t* v,
                           real_t* local_values);

// Copies the locally-stored data in this vector into the given array. 
// The array is assumed to be large enough to store data for all rows of the 
// vector stored locally.
void krylov_vector_copy_out(krylov_vector_t* v,
                            real_t* local_values);

// Assembles added/inserted values into the vector, allowing all the processes
// a consistent representation of the vector. This should be called after calls
// to krylov_vector_set_values/krylov_vector_add_values, and should be 
// placed in between sets and adds. It should also be used after 
// krylov_vector_copy_in.
void krylov_vector_assemble(krylov_vector_t* A);

// Retrieves the values of the elements in the vector identified by the 
// given (global) indices, storing them in the values array. The values 
// must exist on the local process.
void krylov_vector_get_values(krylov_vector_t* v,
                              size_t num_values,
                              index_t* indices,
                              real_t* values);

// Computes and returns the dot product of the vector v with the vector w.
// This is collective, and must be called by all MPI processes.
real_t krylov_vector_dot(krylov_vector_t* v, krylov_vector_t* w);

// Computes and returns the p norm for this vector, where p can be 
// 0 (infinity/max norm), 1, or 2.
// This is collective, and must be called by all MPI processes.
real_t krylov_vector_norm(krylov_vector_t* v, int p);
 
// Computes the weighted 2-norm for this vector, using the weights
// in the vector w. If w is NULL, the 2-norm is returned.
// The weighted 2-norm is sqrt(sum(i, Wi*vi)**2).
// This is collective, and must be called by all MPI processes.
real_t krylov_vector_w2_norm(krylov_vector_t* v, krylov_vector_t* w);

// Computes and returns a weighted root-mean-squared (WRMS) norm for 
// this vector, using the weights in the given vector w. w must be non-NULL.
// The WRMS norm is sqrt(sum(i, (wi*vi)**2)/N).
// This is collective, and must be called by all MPI processes.
real_t krylov_vector_wrms_norm(krylov_vector_t* v, krylov_vector_t* w);
 
// Writes a text representation of the vector (or portion stored on the local
// MPI process) to the given stream.
void krylov_vector_fprintf(krylov_vector_t* v,
                           FILE* stream);

#endif

