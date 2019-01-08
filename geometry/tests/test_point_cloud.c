// Copyright (c) 2012-2019, Jeffrey N. Johnson
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
#include "geometry/point_cloud.h"

static void test_ctor1(void** state)
{
  point_cloud_t* cloud = point_cloud_new(MPI_COMM_SELF, 1000);
  assert_int_equal(1000, cloud->num_points);
  point_cloud_fprintf(cloud, stdout);
  point_cloud_free(cloud);
}

static void test_ctor2(void** state)
{
  int num_points = 1000;
  point_t points[num_points];
  real_t dx = 1.0/num_points;
  for (int i = 0; i < num_points; ++i)
  {
    points[i].x = dx * i;
    points[i].y = 0.0;
    points[i].z = 0.0;
  }
  point_cloud_t* cloud = point_cloud_from_points(MPI_COMM_SELF, points, num_points);
  assert_int_equal(num_points, cloud->num_points);
  point_cloud_fprintf(cloud, stdout);
  point_cloud_free(cloud);
}

static void test_tags(void** state)
{
  point_cloud_t* cloud = point_cloud_new(MPI_COMM_SELF, 1000);
  assert_int_equal(1000, cloud->num_points);

  // Create 10 tags across the index space.
  for (int i = 0; i < 10; ++i)
  {
    char tag_name[129];
    snprintf(tag_name, 128, "tag_%d", i);
    int* tag = point_cloud_create_tag(cloud, tag_name, 100);
    for (int j = 0; j < 100; ++j)
      tag[j] = 100*i + j;
    assert_true(point_cloud_has_tag(cloud, tag_name));
    size_t N;
    int* tag1 = point_cloud_tag(cloud, tag_name, &N);
    assert_true(tag1 == tag);
    assert_int_equal(100, N);
  }

  // Rename the tags, removing the underscore.
  for (int i = 0; i < 10; ++i)
  {
    char old_tag_name[129], new_tag_name[129];
    snprintf(old_tag_name, 128, "tag_%d", i);
    snprintf(new_tag_name, 128, "tag%d", i);
    point_cloud_rename_tag(cloud, old_tag_name, new_tag_name);
    assert_false(point_cloud_has_tag(cloud, old_tag_name));
    assert_true(point_cloud_has_tag(cloud, new_tag_name));
  }

  // Delete the tags.
  for (int i = 0; i < 10; ++i)
  {
    char tag_name[129];
    snprintf(tag_name, 128, "tag%d", i);
    point_cloud_delete_tag(cloud, tag_name);
  }

  point_cloud_free(cloud);
}

static point_cloud_t* cloud_lattice_new(bbox_t* bbox,  
                                        int ni, int nj, int nk)
{
  int N = ni*nj*nk;
  point_cloud_t* cloud = point_cloud_new(MPI_COMM_SELF, N);
  real_t Lx = bbox->x2 - bbox->x1;
  real_t Ly = bbox->y2 - bbox->y1;
  real_t Lz = bbox->z2 - bbox->z1;
  real_t dx = Lx / (ni-1), dy = Ly / (nj-1), dz = Lz / (nk-1);
  int l = 0;
  for (int i = 0; i < ni; ++i)
  {
    real_t xi = bbox->x1 + i*dx;
    for (int j = 0; j < nj; ++j)
    {
      real_t yj = bbox->y1 + j*dy;
      for (int k = 0; k < nk; ++k, ++l)
      {
        real_t zk = bbox->z1 + k*dz;
        cloud->points[l].x = xi;
        cloud->points[l].y = yj;
        cloud->points[l].z = zk;
      }
    }
  }

  return cloud;
}

static void test_unite(void** state)
{
  bbox_t bbox1 = {.x1 = 0.0, .x2 = 1.0,
                  .y1 = 0.0, .y2 = 1.0,
                  .z1 = 0.0, .z2 = 1.0};
  point_cloud_t* c1 = cloud_lattice_new(&bbox1, 10, 10, 10);

  bbox_t bbox2 = {.x1 = 1.0, .x2 = 2.0,
                  .y1 = 0.0, .y2 = 1.0,
                  .z1 = 0.0, .z2 = 1.0};
  point_cloud_t* c2 = cloud_lattice_new(&bbox2, 10, 10, 10);
  point_cloud_unite(c1, c2, "new");
  assert_int_equal(c1->num_points, 2000);
  assert_int_equal(c1->num_ghosts, 2*c2->num_ghosts);
  assert_true(point_cloud_has_tag(c1, "new"));
  size_t N;
  int* new = point_cloud_tag(c1, "new", &N);
  assert_true(new != NULL);
  assert_int_equal(1000, N);

  bbox_t bbox3 = {.x1 = 0.0, .x2 = 1.0,
                  .y1 = 1.0, .y2 = 2.0,
                  .z1 = 0.0, .z2 = 1.0};
  point_cloud_t* c3 = cloud_lattice_new(&bbox3, 10, 10, 10);
  point_cloud_unite(c1, c3, "new");
  assert_int_equal(c1->num_points, 3000);
  assert_int_equal(c1->num_ghosts, 3*c2->num_ghosts);
  assert_true(point_cloud_has_tag(c1, "new"));
  new = point_cloud_tag(c1, "new", &N);
  assert_true(new != NULL);
  assert_int_equal(2000, N);

  point_cloud_free(c1);
  point_cloud_free(c2);
  point_cloud_free(c3);
}

static void test_intersect(void** state)
{
  bbox_t bbox1 = {.x1 = 0.0, .x2 = 1.0,
                  .y1 = 0.0, .y2 = 1.0,
                  .z1 = 0.0, .z2 = 1.0};
  point_cloud_t* c1 = cloud_lattice_new(&bbox1, 11, 11, 11);

  bbox_t bbox2 = {.x1 = 0.3, .x2 = 0.7,
                  .y1 = 0.3, .y2 = 0.7,
                  .z1 = 0.3, .z2 = 0.7};
  point_cloud_t* c2 = cloud_lattice_new(&bbox2, 5, 5, 5);
  point_cloud_intersect(c1, c2, 1e-10);
  assert_int_equal(c1->num_points, 125);

  point_cloud_free(c1);
  point_cloud_free(c2);
}

static void test_difference(void** state)
{
  bbox_t bbox1 = {.x1 = 0.0, .x2 = 1.0,
                  .y1 = 0.0, .y2 = 1.0,
                  .z1 = 0.0, .z2 = 1.0};
  point_cloud_t* c1 = cloud_lattice_new(&bbox1, 11, 11, 11);

  bbox_t bbox2 = {.x1 = 0.3, .x2 = 0.7,
                  .y1 = 0.3, .y2 = 0.7,
                  .z1 = 0.3, .z2 = 0.7};
  point_cloud_t* c2 = cloud_lattice_new(&bbox2, 5, 5, 5);
  point_cloud_difference(c1, c2, 1e-10);
  assert_int_equal(c1->num_points, 11*11*11-5*5*5);

  point_cloud_free(c1);
  point_cloud_free(c2);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_ctor1),
    cmocka_unit_test(test_ctor2),
    cmocka_unit_test(test_tags),
    cmocka_unit_test(test_unite),
    cmocka_unit_test(test_intersect),
    cmocka_unit_test(test_difference)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
