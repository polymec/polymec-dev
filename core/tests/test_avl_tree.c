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
#include "core/avl_tree.h"

int int_values[] = {0, 1, 2, 3, 4, 5};
real_t real_values[] = {0., 1., 2., 3., 4., 5.};

#define DEFINE_AVL_TREE_TEST(tree_name, element) \
void test_##tree_name##_ctor(void** state) \
{ \
  tree_name##_t* t = tree_name##_new(); \
  assert_int_equal(0, tree_name##_size(t)); \
  assert_true(t->root == NULL); \
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
DEFINE_AVL_TREE_TEST(real_avl_tree, real)

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_int_tree),
    unit_test(test_real_tree)
  };
  return run_tests(tests);
}
