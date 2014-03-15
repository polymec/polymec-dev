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
#include "core/mesh.h"

void test_single_cell_mesh_no_topo(void** state)
{
  // Create a single hexahedron without topology.
  mesh_t* mesh = mesh_new(MPI_COMM_WORLD, 1, 0, 6, 8);
  assert_int_equal(1, mesh->num_cells);
  assert_int_equal(0, mesh->num_ghost_cells);
  assert_int_equal(6, mesh->num_faces);
  assert_int_equal(0, mesh->num_edges);
  assert_int_equal(8, mesh->num_nodes);
  mesh_free(mesh);
}

void test_single_cell_mesh_serialization(void** state)
{
  // Create a single hexahedron without topology.
  mesh_t* mesh1 = mesh_new(MPI_COMM_WORLD, 1, 0, 6, 8);
  serializer_t* s = mesh_serializer();
  byte_array_t* buffer = byte_array_new();
  size_t offset = 0;
  serializer_write(s, mesh1, buffer, &offset);
  offset = 0;
  mesh_t* mesh2 = serializer_read(s, buffer, &offset);
  byte_array_free(buffer);
  assert_int_equal(1, mesh2->num_cells);
  assert_int_equal(0, mesh2->num_ghost_cells);
  assert_int_equal(6, mesh2->num_faces);
  assert_int_equal(8, mesh2->num_nodes);
  mesh_free(mesh2);
  mesh_free(mesh1);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_single_cell_mesh_no_topo),
    unit_test(test_single_cell_mesh_serialization),
  };
  return run_tests(tests);
}
