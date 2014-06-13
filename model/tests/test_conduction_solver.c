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
#include "core/silo_file.h"
#include "geometry/create_uniform_mesh.h"
#include "model/conduction_solver.h"

void test_gmres_conduction_solver_ctor(void** state)
{
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  mesh_t* mesh = create_uniform_mesh(MPI_COMM_WORLD, 10, 10, 10, &bbox);
  krylov_solver_t* solver = gmres_conduction_solver_new(mesh, NULL, NULL, 15, 5);
  assert_true(solver != NULL);
  krylov_solver_free(solver);
  mesh_free(mesh);
}

void test_gmres_conduction_solver_laplace_dirichlet(void** state)
{
  // Set up a solver.
  bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  int nx = 10, ny = 10, nz = 10;
  mesh_t* mesh = create_uniform_mesh(MPI_COMM_WORLD, nx, ny, nz, &bbox);
  tag_rectilinear_mesh_faces(mesh, nx, ny, nz, "-x", "+x", "-y", "+y", "-z", "+z");
  krylov_solver_t* solver = gmres_conduction_solver_new(mesh, NULL, NULL, 15, 5);

  // Set Dirichlet boundary conditions.
  const char* tag_names[6] = {"-x", "+x", "-y", "+y", "-z", "+z"};
  real_t alpha[6][nx*ny], gamma[6][nx*ny];
  int num_faces, *tag;
  for (int d = 0; d < 3; ++d)
  {
    tag = mesh_tag(mesh->face_tags, tag_names[2*d], &num_faces);
    for (int i = 0; i < num_faces; ++i)
    {
      alpha[2*d][i] = 1.0;
      gamma[2*d][i] = 0.0;
    }
    conduction_solver_add_bc(solver, tag_names[2*d], alpha[2*d], NULL, gamma[2*d]);

    tag = mesh_tag(mesh->face_tags, tag_names[2*d+1], &num_faces);
    for (int i = 0; i < num_faces; ++i)
    {
      alpha[2*d+1][i] = 1.0;
      gamma[2*d+1][i] = 1.0;
    }
    conduction_solver_add_bc(solver, tag_names[2*d+1], alpha[2*d+1], NULL, gamma[2*d+1]);
  }

  // Find the solution.
  real_t res_norm;
  real_t* X = krylov_solver_vector(solver);
  int num_iters, num_precond;
  bool result = krylov_solver_solve(solver, X, &res_norm, &num_iters, &num_precond);
//  printf("phi = [");
//  for (int i = 0; i < nx*ny*nz; ++i)
//    printf("%g ", X[i]);
//  printf("]\n");
  printf("res_norm = %g, num_iters = %d\n", res_norm, num_iters);
  assert_true(result);

  silo_file_t* silo = silo_file_new(mesh->comm, "test_conduction_solver_laplace_dirichlet", ".", 1, 0, 0, 0.0);
  silo_file_write_mesh(silo, "mesh", mesh);
  silo_file_write_scalar_cell_field(silo, "phi", "mesh", X);
  silo_file_close(silo);

  // Clean up.
  polymec_free(X);
  krylov_solver_free(solver);
  mesh_free(mesh);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_gmres_conduction_solver_ctor),
    unit_test(test_gmres_conduction_solver_laplace_dirichlet)
  };
  return run_tests(tests);
}
