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
#include "geometry/create_tetgen_mesh.h"

void test_create_tetgen_mesh(void** state)
{
  // Create a TetGen mesh from the tetgen_example.* files.
  mesh_t* mesh = create_tetgen_mesh(MPI_COMM_WORLD, 
                                    CMAKE_CURRENT_SOURCE_DIR "/tetgen_example.1.node", 
                                    CMAKE_CURRENT_SOURCE_DIR "/tetgen_example.1.ele", 
                                    CMAKE_CURRENT_SOURCE_DIR "/tetgen_example.1.face", 
                                    CMAKE_CURRENT_SOURCE_DIR "/tetgen_example.1.neigh", 
                                    1);
  mesh_verify(mesh);
  assert_int_equal(1020, mesh->num_cells);
  assert_int_equal(0, mesh->num_ghost_cells);
  assert_int_equal(2286, mesh->num_faces);
  assert_int_equal(1569, mesh->num_edges);
  assert_int_equal(304, mesh->num_nodes);
  mesh_free(mesh);
}

void test_plot_tetgen_mesh(void** state)
{
  // Create a TetGen mesh from the tetgen_example.* files.
  mesh_t* mesh = create_tetgen_mesh(MPI_COMM_WORLD, 
                                    CMAKE_CURRENT_SOURCE_DIR "/tetgen_example.1.node", 
                                    CMAKE_CURRENT_SOURCE_DIR "/tetgen_example.1.ele", 
                                    CMAKE_CURRENT_SOURCE_DIR "/tetgen_example.1.face", 
                                    CMAKE_CURRENT_SOURCE_DIR "/tetgen_example.1.neigh",
                                    1);
  // Plot it.
  silo_file_t* silo = silo_file_new(mesh->comm, "tetgen_example", ".", 1, 0, 0, 0.0);
  silo_file_write_mesh(silo, "mesh", mesh);
  silo_file_close(silo);

  // Clean up.
  mesh_free(mesh);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_create_tetgen_mesh),
    unit_test(test_plot_tetgen_mesh)
  };
  return run_tests(tests);
}
