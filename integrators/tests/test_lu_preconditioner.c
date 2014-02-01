// Copyright (c) 2012-2014, Jeffrey N. Johnson
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this 
// list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice, 
// this list of conditions and the following disclaimer in the documentation 
// and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>

#include "cmockery.h"
#include "slu_ddefs.h"
#include "core/polymec.h"
#include "core/array_utils.h"
#include "geometry/create_uniform_mesh.h"
#include "integrators/lu_preconditioners.h"

static adj_graph_t* graph_from_uniform_mesh()
{
  bbox_t box = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  mesh_t* m = create_uniform_mesh(MPI_COMM_WORLD, 10, 10, 10, &box);
  adj_graph_t* g = graph_from_mesh_cells(m);
  mesh_free(m);
  return g;
}

static int sys_func(void* context, real_t t, real_t* x, real_t* F)
{
  return 0;
}

void test_ctor(void** state)
{
  adj_graph_t* g = graph_from_uniform_mesh();
  preconditioner_t* precond = lu_preconditioner_new(NULL, sys_func, NULL, g);
  preconditioner_free(precond);
  adj_graph_free(g);
}

void test_matrix(void** state)
{
  // Build a matrix A.
  adj_graph_t* g = graph_from_uniform_mesh();
  preconditioner_t* precond = lu_preconditioner_new(NULL, sys_func, NULL, g);
  preconditioner_matrix_t* mat = preconditioner_matrix(precond);

  // Check the non-zero structure of A.
  SuperMatrix* A = preconditioner_matrix_context(mat);
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
  preconditioner_matrix_free(mat);
  preconditioner_free(precond);
  adj_graph_free(g);
}

/*******************************************************************************
**
** Example from Dennis and Schnabel (1996) Numerical Methods for
** Unconstrained Optimization and Nonlinear Equations, pg 87
**
*******************************************************************************/
typedef struct {
  int num_unknown;
} context_t;

static int dennis_schnabel_1(void* context, real_t t, real_t* x, real_t* F) {
  // F(x) = [ x1 + x2 - 3.0,
  //        [ x1^2 + x2^2 - 9.0

  F[0] = x[0] + x[1] - 3.0;
  F[1] = x[0] * x[0] + x[1] * x[1] - 9.0;
  return 0;
}

void test_numerical_jacobian_ds1(void **state) {
  adj_graph_t* g = adj_graph_new(MPI_COMM_WORLD, 1);
  adj_graph_t* bg = adj_graph_new_with_block_size(2, g);
  adj_graph_free(g);
  context_t context;
  preconditioner_t* precond = lu_preconditioner_new(&context, dennis_schnabel_1, NULL, bg);

  real_t time = 0.0;
  real_t x[2];
  x[0] = 1.0;
  x[1] = 5.0;
  preconditioner_matrix_t* mat = preconditioner_matrix(precond);
  preconditioner_compute_jacobian(precond, time, x, mat);

  // expected jacobian:
  // J = [ 1  1  ]
  //     [ 2  10 ]
  assert_true(fabs(preconditioner_matrix_coeff(mat, 0, 0) - 1.0) < 1e-14);
  assert_true(fabs(preconditioner_matrix_coeff(mat, 0, 1) - 1.0) < 1e-14);
  assert_true(fabs(preconditioner_matrix_coeff(mat, 1, 0) - 2.0) < 1e-14);
  assert_true(fabs(preconditioner_matrix_coeff(mat, 1, 1) - 10.0) < 1e-14);
  preconditioner_matrix_free(mat);
  adj_graph_free(bg);
  preconditioner_free(precond);
}  // end test_numerical_jacobian

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_ctor),
    unit_test(test_matrix),
    unit_test(test_numerical_jacobian_ds1)
  };
  return run_tests(tests);
}
