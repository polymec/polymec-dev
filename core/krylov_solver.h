// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_KRYLOV_SOLVER_H
#define POLYMEC_KRYLOV_SOLVER_H

#include "core/polymec.h"
#include "core/adj_graph.h"
#include "core/unordered_map.h"

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

// Objects of this type construct Krylov matrices, vectors, and solvers. 
// A Krylov factory exposes capabilities within a given library (PETSc, HYPRE, 
// etc). The factory must continue to exist as long as any of the solvers, 
// matrices, and vectors it has created.
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
  krylov_matrix_t* (*matrix)(void* context, adj_graph_t* sparsity);
  krylov_matrix_t* (*block_matrix)(void* context, adj_graph_t* sparsity, int block_size);
  krylov_vector_t* (*vector)(void* context, adj_graph_t* dist_graph);
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
  void* (*clone)(void* context);
  void (*zero)(void* context);
  void (*scale)(void* context, real_t scale_factor);
  void (*add_identity)(void* context, real_t scale_factor);
  void (*add_diagonal)(void* context, void* D);
  void (*set_diagonal)(void* context, void* D);
  void (*set_values)(void* context, index_t num_rows, index_t* num_columns, index_t* rows, index_t* columns, real_t* values);
  void (*add_values)(void* context, index_t num_rows, index_t* num_columns, index_t* rows, index_t* columns, real_t* values);
  void (*start_assembly)(void* context);
  void (*finish_assembly)(void* context);
  void (*get_values)(void* context, index_t num_rows, index_t* num_columns, index_t* rows, index_t* columns, real_t* values);
  void (*dtor)(void* context);
} krylov_matrix_vtable;

// This constructor should be called with a context pointer and a virtual 
// table to create an instance of a krylov_matrix subclass.
krylov_matrix_t* krylov_matrix_new(void* context,
                                   krylov_matrix_vtable vtable,
                                   int num_local_rows,
                                   index_t num_global_rows);

// This virtual table must be filled out for any subclass of krylov_matrix.
typedef struct
{
  void* (*clone)(void* context);
  void (*zero)(void* context);
  void (*set_value)(void* context, real_t value);
  void (*scale)(void* context, real_t scale_factor);
  void (*set_values)(void* context, index_t num_values, index_t* indices, real_t* values);
  void (*add_values)(void* context, index_t num_values, index_t* indices, real_t* values);
  void (*start_assembly)(void* context);
  void (*finish_assembly)(void* context);
  void (*get_values)(void* context, index_t num_values, index_t* indices, real_t* values);
  real_t (*norm)(void* context, int p);
  void (*dtor)(void* context);
} krylov_vector_vtable;

// This constructor should be called with a context pointer and a virtual 
// table to create an instance of a krylov_vector subclass.
krylov_vector_t* krylov_vector_new(void* context,
                                   krylov_vector_vtable vtable,
                                   int local_size,
                                   index_t global_size);

//------------------------------------------------------------------------
//                  Bundled third-party Krylov factories 
//------------------------------------------------------------------------

// This creates a PETSc-based Krylov factory that can be used for constructing
// matrices, vectors, solvers, using the given petsc directory and architecture 
// (or the environment variables PETSC_DIR and PETSC_ARCH if these strings are 
// NULL) to find the underlying PETSc implementation. If no such 
// underlying implementation can be found, this function returns NULL.
krylov_factory_t* petsc_krylov_factory(const char* petsc_dir,
                                       const char* petsc_arch);

// This creates a HYPRE-based Krylov factory that can be used for constructing
// matrices, vectors, solvers, using the HYPRE library (libHYPRE.so, etc) 
// located in the given directory. If no such library can be found or loaded, 
// this function returns NULL.
krylov_factory_t* hypre_krylov_factory(const char* hypre_dir);

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
                                       adj_graph_t* sparsity);

// Constructs a (square) Krylov sparse block matrix with the sparsity pattern 
// given by an adjacency graph, and the given block size.
krylov_matrix_t* krylov_factory_block_matrix(krylov_factory_t* factory, 
                                             adj_graph_t* sparsity,
                                             int block_size);

// Constructs a vector on the given communicator with its local and global 
// dimensions defined by the given distributed graph. This graph is typically 
// the same as the sparsity graph passed to a corresponding matrix.
krylov_vector_t* krylov_factory_vector(krylov_factory_t* factory,
                                       adj_graph_t* dist_graph);

// Constructs a preconditioned conjugate gradient (PCG) Krylov solver. Keep in 
// mind that this method can only be used for systems having symmetric, 
// positive-definite matrices.
krylov_solver_t* krylov_factory_pcg_solver(krylov_factory_t* factory,
                                           MPI_Comm comm);
                                             
// Constructs a GMRES Krylov solver with the given Krylov subspace dimension.
krylov_solver_t* krylov_factory_gmres_solver(krylov_factory_t* factory,
                                             MPI_Comm comm,
                                             int krylov_dimension);
                                             
// Constructs a stabilized bi-conjugate Krylov solver.
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
// in the vector b. Returns true if the solution was obtained, false if not.
// The operator matrix A must be set using krylov_solver_set_operator before
// this function is called. The number of linear iterations will be stored in 
// num_iterations upon success. If the solve is unsuccessful, the vector b 
// will NOT contain the solution unless the residual was reduced.
bool krylov_solver_solve(krylov_solver_t* solver, 
                         krylov_vector_t* b, 
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

// Creates and returns a deep copy of a matrix.
krylov_matrix_t* krylov_matrix_clone(krylov_matrix_t* A);

// Returns a pointer to the underlying matrix implementation. You can use 
// this if you have explicitly linked your program to the library providing
// the implementation.
void* krylov_matrix_impl(krylov_matrix_t* A);

// Returns the number of locally stored rows in the matrix.
int krylov_matrix_num_local_rows(krylov_matrix_t* A);

// Returns the number of globally stored rows in the matrix.
int krylov_matrix_num_global_rows(krylov_matrix_t* A);

// Zeros all of the entries in the given matrix.
void krylov_matrix_zero(krylov_matrix_t* A);

// Multiplies this matrix by a scale factor.
void krylov_matrix_scale(krylov_matrix_t* A,
                         real_t scale_factor);

// Adds a scaled identity matrix to this one.
void krylov_matrix_add_identity(krylov_matrix_t* A,
                                real_t scale_factor);

// Adds the entries of the given vector (D) to the diagonal of the given matrix (A).
void krylov_matrix_add_diagonal(krylov_matrix_t* A,
                                krylov_vector_t* D);

// Sets the diagonal entries in the given matrix (A) to those in the given vector (D).
void krylov_matrix_set_diagonal(krylov_matrix_t* A,
                                krylov_vector_t* D);

// Sets the values of the elements in the matrix identified by the given 
// (globally-indexed) rows and columns. The rows and columns being set must 
// exist on the local process. num_columns and rows are arrays of size num_rows, 
// and contain the numbers of columns and the row indices for each row. columns
// is an array containing the indices of all the columns, in row-major, column-minor
// order. values is an array with values corresponding to the entries in the columns 
// array.
void krylov_matrix_set_values(krylov_matrix_t* A,
                              index_t num_rows,
                              index_t* num_columns,
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
                              index_t num_rows,
                              index_t* num_columns,
                              index_t* rows, index_t* columns,
                              real_t* values);
                              
// Assemble the matrix. This should be called after setting or adding to its 
// values.
void krylov_matrix_assemble(krylov_matrix_t* A);

// Start the (possibly asynchronous) matrix assembly. This should be called 
// after setting or adding to its values.
void krylov_matrix_start_assembly(krylov_matrix_t* A);

// Finish the (possibly asynchronous) matrix assembly. This should be called 
// after setting or adding to its values.
void krylov_matrix_finish_assembly(krylov_matrix_t* A);

// Retrieves the values of the elements in the matrix identified by the 
// given (globally-indexed) rows and columns, storing them in the values array.
void krylov_matrix_get_values(krylov_matrix_t* A,
                              index_t num_rows,
                              index_t* num_columns,
                              index_t* rows, index_t* columns,
                              real_t* values);

//------------------------------------------------------------------------
//                      Krylov vector interface
//------------------------------------------------------------------------

// Frees a vector.
void krylov_vector_free(krylov_vector_t* v);

// Creates and returns a deep copy of a vector.
krylov_vector_t* krylov_vector_clone(krylov_vector_t* v);

// Returns a pointer to the underlying vector implementation. You can use 
// this if you have explicitly linked your program to the library providing
// the implementation.
void* krylov_vector_impl(krylov_vector_t* v);

// Returns the locally-stored size (dimension) of the vector.
int krylov_vector_local_size(krylov_vector_t* v);

// Returns the global size (dimension) of the vector.
int krylov_vector_global_size(krylov_vector_t* v);

// Zeros all of the entries in the given vector.
void krylov_vector_zero(krylov_vector_t* v);

// Sets all of the entries in the given vector to the given value.
void krylov_vector_set_value(krylov_vector_t* v,
                             real_t value);

// Scales the vector by the given factor.
void krylov_vector_scale(krylov_vector_t* v,
                         real_t scale_factor);

// Sets the values of the elements in the vector identified by the given 
// (globally-indexed) indices.
void krylov_vector_set_values(krylov_vector_t* v,
                              index_t num_values,
                              index_t* indices,
                              real_t* values);
                              
// Adds the given values of the elements to those in the vector.
// The values are identified by the given (globally-indexed) indices.
void krylov_vector_add_values(krylov_vector_t* v,
                              index_t num_values,
                              index_t* indices,
                              real_t* values);
                              
// Assemble the vector. This should be called after setting or adding to its 
// values.
void krylov_vector_assemble(krylov_vector_t* v);

// Starts the (possibly asynchronous) vector assembly. This should be called 
// after setting or adding to its values.
void krylov_vector_start_assembly(krylov_vector_t* v);

// Finish the (possibly asynchronous) matrix assembly. This should be called 
// after setting or adding to its values.
void krylov_vector_finish_assembly(krylov_vector_t* v);

// Retrieves the values of the elements in the vector identified by the 
// given (global) indices, storing them in the values array. The values 
// must exist on the local process.
void krylov_vector_get_values(krylov_vector_t* v,
                              index_t num_values,
                              index_t* indices,
                              real_t* values);

// Computes and returns the Lp norm for this vector, where p can be 
// 0 (infinity/max norm), 1, or 2.
real_t krylov_vector_norm(krylov_vector_t* v, int p);
 
#endif

