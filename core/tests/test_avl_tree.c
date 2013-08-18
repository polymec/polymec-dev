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
#include "core/avl_tree.h"

int int_values[] = {0, 1, 2, 3, 4, 5};
double double_values[] = {0., 1., 2., 3., 4., 5.};

#define DEFINE_AVL_TREE_TEST(tree_name, element) \
void test_##tree_name##_ctor(void** state) \
{ \
  tree_name##_t* t = tree_name##_new(); \
  assert_int_equal(0, tree_name##_size(t)); \
  assert_true(t->root == NULL); \
  assert_true(t->arena == NULL); \
  tree_name##_free(t); \
} \
\
void test_##tree_name##_insert(void** state) \
{ \
  tree_name##_t* t = tree_name##_new(); \
  for (int i = 0; i < 5; ++i) \
    tree_name##_insert(t, element##_values[i]); \
  assert_int_equal(5, tree_name##_size(t)); \
  assert_true(t->root != NULL); \
  tree_name##_clear(t); \
  assert_int_equal(0, tree_name##_size(t)); \
  assert_true(t->root == NULL); \
  tree_name##_free(t); \
} \
\
void test_##element##_tree(void** state) \
{ \
  test_##tree_name##_ctor(state); \
  test_##tree_name##_insert(state); \
}

DEFINE_AVL_TREE_TEST(int_avl_tree, int)
DEFINE_AVL_TREE_TEST(double_avl_tree, double)

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_int_tree),
    unit_test(test_double_tree)
  };
  return run_tests(tests);
}
