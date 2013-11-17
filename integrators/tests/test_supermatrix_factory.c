// Copyright 2012-2013 Jeffrey Johnson.
// 
// This file is part of Polymec, and is licensed under the Apache License, 
// Version 2.0 (the "License"); you may not use this file except in 
// compliance with the License. You may may find the text of the license in 
// the LICENSE file at the top-level source directory, or obtain a copy of 
// it at
// 
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmockery.h"
#include "core/polymec.h"
#include "core/graph_from_mesh_cells.h"
#include "core/create_uniform_mesh.h"
#include "core/array_utils.h"
#include "integrators/supermatrix_factory.h"

static adj_graph_t* graph_from_uniform_mesh()
{
  int nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  bbox_t box = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  mesh_t* m = create_uniform_mesh(MPI_COMM_WORLD, 10, 10, 10, &box);
  adj_graph_t* g = graph_from_mesh_cells(m);
  mesh_free(m);
  return g;
}

static int sys_func(N_Vector u, N_Vector F, void* context)
{
  return 0;
}

typedef struct
{
  double time;
} sys_func_t;

static void set_sys_func_time(double t, void* sys_func)
{
  sys_func_t* F = sys_func;
  F->time = t;
}

void test_F_ctor(void** state)
{
  adj_graph_t* g = graph_from_uniform_mesh();

  sys_func_t F = {.time = 0.0};
  supermatrix_factory_t* factory = 
    supermatrix_factory_from_sys_func(g, sys_func, NULL, &F);

  supermatrix_factory_free(factory);
  adj_graph_free(g);
}

void test_time_dep_F_ctor(void** state)
{
  adj_graph_t* g = graph_from_uniform_mesh();

  sys_func_t F = {.time = 0.0};
  supermatrix_factory_t* factory = 
    supermatrix_factory_from_sys_func(g, sys_func, set_sys_func_time, &F);

  supermatrix_factory_free(factory);
  adj_graph_free(g);
}

static int rhs_func(double t, N_Vector u, N_Vector F, void* context)
{
  return 0;
}

typedef struct
{
} rhs_t;

void test_rhs_ctor(void** state)
{
  adj_graph_t* g = graph_from_uniform_mesh();

  rhs_t rhs = {};
  supermatrix_factory_t* factory = 
    supermatrix_factory_from_rhs(g, rhs_func, &rhs);

  supermatrix_factory_free(factory);
  adj_graph_free(g);
}

void test_supermatrix(void** state)
{
  // Build a supermatrix A.
  adj_graph_t* g = graph_from_uniform_mesh();
  sys_func_t F = {.time = 0.0};
  supermatrix_factory_t* factory = 
    supermatrix_factory_from_sys_func(g, sys_func, NULL, &F);
  SuperMatrix* A = supermatrix_factory_matrix(factory);

  // Check the non-zero structure of A.
  int N = adj_graph_num_vertices(g);
  NRformat* Adata = A->Store;
  for (int i = 0; i < N; ++i)
  {
    int row_index = Adata->rowptr[i];
    int next_row_index = Adata->rowptr[i+1];
    int num_cols = next_row_index - row_index;

    // The column indices in row i start with the diagonal element
    // and are followed by the off-diagonals in ascending order.

    // Check the diagonal element.
    {
      int j = Adata->colind[row_index];
      assert_int_equal(j, i);
    }
    
    // Check off-diagonals.
    int pos = 0, j;
    while (adj_graph_next_edge(g, i, &pos, &j))
    {
      int* jptr = int_bsearch(&Adata->colind[row_index+1], num_cols-1, j);
      assert_true(jptr != NULL);
    }
  }

  // Clean up.
  Destroy_CompRow_Matrix(A);
  supermatrix_factory_free(factory);
  adj_graph_free(g);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_F_ctor),
    unit_test(test_time_dep_F_ctor),
    unit_test(test_rhs_ctor),
    unit_test(test_supermatrix)
  };
  return run_tests(tests);
}
