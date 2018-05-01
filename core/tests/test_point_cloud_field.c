// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmocka.h"
#include "core/point_cloud_field.h"

static void test_ctor(void** state)
{
  point_cloud_t* cloud = point_cloud_new(MPI_COMM_SELF, 1000);
  point_cloud_field_t* f1 = point_cloud_field_new(cloud, 1);
  point_cloud_field_t* f2 = point_cloud_field_new(cloud, 2);
  point_cloud_field_t* f3 = point_cloud_field_new(cloud, 3);
  assert_int_equal(point_cloud_field_num_components(f1), 1);
  assert_int_equal(point_cloud_field_num_local_values(f1), 1000);
  assert_int_equal(point_cloud_field_num_ghost_values(f1), 0);
  assert_int_equal(point_cloud_field_num_components(f2), 2);
  assert_int_equal(point_cloud_field_num_local_values(f2), 1000);
  assert_int_equal(point_cloud_field_num_ghost_values(f2), 0);
  assert_int_equal(point_cloud_field_num_components(f3), 3);
  assert_int_equal(point_cloud_field_num_local_values(f3), 1000);
  assert_int_equal(point_cloud_field_num_ghost_values(f3), 0);
  point_cloud_field_free(f1);
  point_cloud_field_free(f2);
  point_cloud_field_free(f3);
  point_cloud_free(cloud);
}

static void test_set_ghosts(void** state)
{
  point_cloud_t* cloud = point_cloud_new(MPI_COMM_SELF, 1000);
  point_cloud_field_t* f = point_cloud_field_new(cloud, 1);
  assert_int_equal(point_cloud_field_num_local_values(f), 1000);
  assert_int_equal(point_cloud_field_num_ghost_values(f), 0);
  point_cloud_set_num_ghosts(cloud, 100);
  assert_int_equal(point_cloud_field_num_local_values(f), 1000);
  assert_int_equal(point_cloud_field_num_ghost_values(f), 100);
  point_cloud_field_free(f);
  point_cloud_free(cloud);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_ctor),
    cmocka_unit_test(test_set_ghosts)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
