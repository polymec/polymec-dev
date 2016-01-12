// Copyright (c) 2012-2016, Jeffrey N. Johnson
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
#include "core/silo_file.h"
#include "geometry/create_uniform_mesh.h"
#include "geometry/create_welded_block_mesh.h"

void test_2_block_weld(void** state)
{
  // Weld two blocks together!
  bbox_t bbox1 = {.x1 = -0.5, .x2 = 0.0,
                  .y1 = -0.5, .y2 = 0.5,
                  .z1 = -0.5, .z2 = 0.5};
  mesh_t* block1 = create_uniform_mesh(MPI_COMM_SELF, 5, 10, 10, &bbox1);
  tag_rectilinear_mesh_faces(block1, 
                             "left", "middle",
                             "front_left", "back_left",
                             "bottom_left", "top_left");

  bbox_t bbox2 = {.x1 =  0.0, .x2 = 0.5,
                  .y1 = -0.5, .y2 = 0.5,
                  .z1 = -0.5, .z2 = 0.5};
  mesh_t* block2 = create_uniform_mesh(MPI_COMM_SELF, 5, 10, 10, &bbox2);
  tag_rectilinear_mesh_faces(block2, 
                             "middle", "right",
                             "front_right", "back_right",
                             "bottom_right", "top_right");

  mesh_t* blocks[2] = {block1, block2};
  mesh_t* weld = create_welded_block_mesh(blocks, 2, 1e-8);
  assert_true(mesh_verify_topology(weld, polymec_error));
  assert_true(weld->comm == MPI_COMM_SELF);
  assert_int_equal(1000, weld->num_cells);
  assert_int_equal(3300, weld->num_faces);
  assert_int_equal(1331, weld->num_nodes);

  // Plot it!
  silo_file_t* silo = silo_file_new(MPI_COMM_SELF, "2_block_weld", "", 1, 0, 0, 0.0);
  silo_file_write_mesh(silo, "mesh", weld);
  silo_file_close(silo);

  // Clean up.
  mesh_free(block1);
  mesh_free(block2);
  mesh_free(weld);
}

void test_L_weld(void** state)
{
  // Weld three blocks together in an L shape!
  bbox_t bbox1 = {.x1 = -0.5, .x2 = 0.0,
                  .y1 = -0.5, .y2 = 0.0,
                  .z1 = -0.25, .z2 = 0.25};
  mesh_t* block1 = create_uniform_mesh(MPI_COMM_SELF, 10, 10, 10, &bbox1);
  tag_rectilinear_mesh_faces(block1, 
                             "left", "middle",
                             "front_left", "middle_left",
                             "bottom_left", "top_left");

  bbox_t bbox2 = {.x1 =  0.0, .x2 = 0.5,
                  .y1 = -0.5, .y2 = 0.0,
                  .z1 = -0.25, .z2 = 0.25};
  mesh_t* block2 = create_uniform_mesh(MPI_COMM_SELF, 10, 10, 10, &bbox2);
  tag_rectilinear_mesh_faces(block2, 
                             "middle", "right",
                             "front_right", "back_right",
                             "bottom_right", "top_right");

  bbox_t bbox3 = {.x1 = -0.5, .x2 = 0.0,
                  .y1 =  0.0, .y2 = 0.5,
                  .z1 = -0.25, .z2 = 0.25};
  mesh_t* block3 = create_uniform_mesh(MPI_COMM_SELF, 10, 10, 10, &bbox3);
  tag_rectilinear_mesh_faces(block3, 
                             "far_left", "far_middle",
                             "middle_left", "back_left",
                             "far_bottom_left", "far_top_left");

  mesh_t* blocks[3] = {block1, block2, block3};
  mesh_t* weld = create_welded_block_mesh(blocks, 3, 1e-8);
  assert_true(mesh_verify_topology(weld, polymec_error));
  assert_int_equal(3000, weld->num_cells);
  assert_true(weld->comm == MPI_COMM_SELF);

  // Plot it!
  silo_file_t* silo = silo_file_new(MPI_COMM_SELF, "L_weld", "", 1, 0, 0, 0.0);
  silo_file_write_mesh(silo, "mesh", weld);
  silo_file_close(silo);

  // Clean up.
  mesh_free(block1);
  mesh_free(block2);
  mesh_free(block3);
  mesh_free(weld);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_2_block_weld),
    unit_test(test_L_weld)
  };
  return run_tests(tests);
}
