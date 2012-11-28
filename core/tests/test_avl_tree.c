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
