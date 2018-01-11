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
#include <time.h>
#include "cmocka.h"
#include "core/octree.h"

static octree_t* tree_with_N_points(int N)
{
  rng_t* rng = host_rng_new();
  bbox_t bounding_box = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  octree_t* tree = octree_new(&bounding_box);
  for (int i = 0; i < N; ++i)
  {
    point_t x;
    point_randomize(&x, rng, &bounding_box);
    octree_insert(tree, &x, i);
  }
  return tree;
}

static void test_insert(void** state) 
{ 
  // Create a tree containing 100 random points.
  octree_t* tree = tree_with_N_points(100);
  assert_int_equal(100, octree_size(tree));
  octree_free(tree);
}

static void test_delete(void** state) 
{ 
  // Create a tree containing 100 random points.
  octree_t* tree = tree_with_N_points(100);

  // Now delete each of them individually.
  for (int i = 0; i < 100; ++i) 
  {
    octree_delete(tree, i);
    assert_int_equal(100-i-1, octree_size(tree));
  }

  octree_free(tree);
}

typedef struct 
{
  int* branch_data;
  int* leaf_data;
  bool passed;
} tree_data_t;

static bool verify_branch_marked(void* context, int depth, int branch_index, int parent_index)
{
  tree_data_t* data = context;

  if (data->branch_data[branch_index] == 0)
    data->passed = false;

  return true; // Visit the kids.
}

static void verify_leaf_marked(void* context, int leaf_index, point_t* point, int parent_index)
{
  tree_data_t* data = context;
  if (data->leaf_data[leaf_index] == 0)
    data->passed = false;
}

static bool mark_branch_pre(void* context, int depth, int branch_index, int parent_index)
{
  tree_data_t* data = context;

  // The parent shouldn't be marked yet.
  ASSERT((parent_index == -1) || (data->branch_data[parent_index] == 0));
  data->branch_data[branch_index] = 1;

  return true; // Visit the kids.
}

static void mark_leaf_pre(void* context, int leaf_index, point_t* point, int parent_index)
{
  tree_data_t* data = context;

  // The parent shouldn't be marked yet.
  ASSERT(data->branch_data[parent_index] == 0);
  data->leaf_data[leaf_index] = 1;
}

static void test_visit_pre(void** state) 
{ 
  // Create a tree containing 100 random points.
  octree_t* tree = tree_with_N_points(100);

  // Create data for the tree.
  tree_data_t data;
  data.branch_data = octree_new_branch_array(tree, sizeof(int));
  data.leaf_data = octree_new_leaf_array(tree, sizeof(int));
  data.passed = true;

  // Visit the nodes in the tree, marking them as visited.
  octree_visit(tree, OCTREE_PRE, &data, mark_branch_pre, mark_leaf_pre);

  // Make sure everything is marked.
  octree_visit(tree, OCTREE_PRE, &data, verify_branch_marked, verify_leaf_marked);
  assert_true(data.passed);

  // Clean up.
  polymec_free(data.leaf_data);
  polymec_free(data.branch_data);
  octree_free(tree);
}

static bool mark_branch_post(void* context, int depth, int branch_index, int parent_index)
{
  tree_data_t* data = context;

  // The parent should already be marked.
  ASSERT((parent_index == -1) || (data->branch_data[parent_index] == 1));
  data->branch_data[branch_index] = 1;

  return true; // Visit the kids.
}

static void mark_leaf_post(void* context, int leaf_index, point_t* point, int parent_index)
{
  tree_data_t* data = context;

  // The parent should already be marked.
  ASSERT(data->branch_data[parent_index] == 1);
  data->leaf_data[leaf_index] = 1;
}

static void test_visit_post(void** state) 
{ 
  // Create a tree containing 100 random points.
  octree_t* tree = tree_with_N_points(100);

  // Create data for the tree.
  tree_data_t data;
  data.branch_data = octree_new_branch_array(tree, sizeof(int));
  data.leaf_data = octree_new_leaf_array(tree, sizeof(int));
  data.passed = true;

  // Visit the nodes in the tree, marking them as visited.
  octree_visit(tree, OCTREE_POST, &data, mark_branch_post, mark_leaf_post);

  // Make sure everything is marked.
  octree_visit(tree, OCTREE_POST, &data, verify_branch_marked, verify_leaf_marked);
  assert_true(data.passed);

  // Clean up.
  polymec_free(data.leaf_data);
  polymec_free(data.branch_data);
  octree_free(tree);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const struct CMUnitTest tests[] = 
  {
    cmocka_unit_test(test_insert),
    cmocka_unit_test(test_delete),
    cmocka_unit_test(test_visit_pre),
    cmocka_unit_test(test_visit_post)
  };
  return cmocka_run_group_tests(tests, NULL, NULL);
}
