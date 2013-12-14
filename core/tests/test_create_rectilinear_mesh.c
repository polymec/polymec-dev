// Copyright (c) 2012-2013, Jeffrey N. Johnson
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
#include "core/write_silo.h"
#include "core/create_rectilinear_mesh.h"
#include "polytope_c.h"

void test_create_rectilinear_mesh(void** state)
{
  // Create a 10x10x10 rectilinear mesh.
  double xs[] = {0.0, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0, 256.0, 512.0};
  double ys[] = {0.0, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0, 256.0, 512.0};
  double zs[] = {0.0, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0, 256.0, 512.0};
  mesh_t* mesh = create_rectilinear_mesh(MPI_COMM_WORLD, xs, 11, ys, 11, zs, 11);
  mesh_verify(mesh);
  assert_int_equal(10*10*10, mesh->num_cells);
  assert_int_equal(0, mesh->num_ghost_cells);
  assert_int_equal(mesh->num_faces, 3*10*10*11);
  assert_int_equal(mesh->num_edges, 3*10*11*11);
  assert_int_equal(11*11*11, mesh->num_nodes);
  mesh_free(mesh);
}

void test_plot_rectilinear_mesh(void** state)
{
  // Create a 4x4x4 rectilinear mesh.
  double xs[] = {0.0, 1.0, 2.0, 4.0, 8.0};
  double ys[] = {0.0, 1.0, 2.0, 4.0, 8.0};
  double zs[] = {0.0, 1.0, 2.0, 4.0, 8.0};
  mesh_t* mesh = create_rectilinear_mesh(MPI_COMM_WORLD, xs, 5, ys, 5, zs, 5);

  // Plot it.
  double ones[4*4*4];
  for (int c = 0; c < 4*4*4; ++c)
    ones[c] = 1.0*c;
  string_ptr_unordered_map_t* fields = string_ptr_unordered_map_new();
  string_ptr_unordered_map_insert(fields, "solution", ones);
  write_silo_mesh(mesh, NULL, NULL, NULL, fields, "rectilinear_4x4x4", ".",
                  0, 0.0, MPI_COMM_SELF, 1, 0);

  // Clean up.
  string_ptr_unordered_map_free(fields);
  mesh_free(mesh);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_create_rectilinear_mesh),
    unit_test(test_plot_rectilinear_mesh)
  };
  return run_tests(tests);
}
