// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>

#include "cmockery.h"
#include "core/polymec.h"
#include "core/file_utils.h"
#include "core/krylov_solver.h"

// This helper creates a distributed graph for a 1D finite difference 
// discretization of the Laplacian operator.
static adj_graph_t* create_1d_laplacian_graph(int N_local)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  int nprocs, rank;
  MPI_Comm_size(comm, &nprocs);
  MPI_Comm_rank(comm, &rank);

  adj_graph_t* graph = adj_graph_new(MPI_COMM_WORLD, N_local);

  // Set interior edges.
  for (int i = 1; i < N_local-1; ++i)
  {
    adj_graph_set_num_edges(graph, i, 2);
    int* edges = adj_graph_edges(graph, i);
    edges[0] = i-1;
    edges[1] = i+1;
  }

  // Set boundary edges.
  if (nprocs == 1)
  {
    adj_graph_set_num_edges(graph, 0, 1);
    int* edges = adj_graph_edges(graph, 0);
    edges[0] = 1;

    adj_graph_set_num_edges(graph, N_local-1, 1);
    edges = adj_graph_edges(graph, N_local-1);
    edges[0] = N_local-2;
  }
  else if (rank == 0)
  {
    adj_graph_set_num_edges(graph, 0, 1);
    int* edges = adj_graph_edges(graph, 0);
    edges[0] = 1;

    adj_graph_set_num_edges(graph, N_local-1, 2);
    edges = adj_graph_edges(graph, N_local-1);
    edges[0] = N_local-2;
    edges[1] = N_local+1;
  }
  else if (rank == (nprocs-1))
  {
    adj_graph_set_num_edges(graph, 0, 2);
    int* edges = adj_graph_edges(graph, 0);
    edges[0] = N_local;
    edges[1] = 1;

    adj_graph_set_num_edges(graph, N_local-1, 1);
    edges = adj_graph_edges(graph, N_local-1);
    edges[0] = N_local-2;
  }
  else
  {
    adj_graph_set_num_edges(graph, 0, 2);
    int* edges = adj_graph_edges(graph, 0);
    edges[0] = N_local;
    edges[1] = 1;

    adj_graph_set_num_edges(graph, N_local-1, 2);
    edges = adj_graph_edges(graph, N_local-1);
    edges[0] = N_local-2;
    edges[1] = N_local+1;
  }

  return graph;
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
    // Create a distributed graph with 1000 local vertices.
    adj_graph_t* graph = create_1d_laplacian_graph(1000);

    // Create a matrix for this graph.
    krylov_matrix_t* mat = krylov_factory_matrix(factory, graph);
    assert_int_equal(1000, krylov_matrix_num_local_rows(mat));
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    assert_int_equal(1000*nprocs, krylov_matrix_num_global_rows(mat));

    // Clone it.
    krylov_matrix_t* mat1 = krylov_matrix_clone(mat);
    assert_int_equal(1000, krylov_matrix_num_local_rows(mat1));

    // Put everything away.
    krylov_matrix_free(mat);
    krylov_matrix_free(mat1);
    krylov_factory_free(factory);
  }
}

static void test_krylov_vector(void** state, krylov_factory_t* factory)
{
  if (factory != NULL)
  {
    // Create a distributed graph with 1000 local vertices.
    adj_graph_t* graph = create_1d_laplacian_graph(1000);

    // Create a vector for this graph.
    krylov_vector_t* vec = krylov_factory_vector(factory, graph);
    assert_int_equal(1000, krylov_vector_local_size(vec));
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
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

void test_petsc_krylov_vector(void** state)
{
  krylov_factory_t* petsc = create_petsc_krylov_factory();
  test_krylov_vector(state, petsc);
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

void test_hypre_krylov_vector(void** state)
{
  krylov_factory_t* hypre = create_hypre_krylov_factory();
  test_krylov_vector(state, hypre);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_petsc_krylov_factory),
    unit_test(test_petsc_krylov_matrix),
    unit_test(test_petsc_krylov_vector),
    unit_test(test_hypre_krylov_factory),
    unit_test(test_hypre_krylov_matrix),
    unit_test(test_hypre_krylov_vector)
  };
  return run_tests(tests);
}

