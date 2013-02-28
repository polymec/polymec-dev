#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmockery.h"
#include "geometry/block.h"
#include "geometry/ball_mapping.h"

static void test_ball_mapping(block_t* block, int block_index)
{
  point_t zero = {.x = 0.0, .y = 0.0, .z = 0.0};
  sp_func_t* mapping = ball_mapping_new(&zero, 1.0, 0.5, block_index);
  int num_points = block_num_points(block);
  point_t xis[num_points], xs[num_points];
  double xjs[3*num_points];
  block_get_points(block, xis);
  for (int i = 0; i < num_points; ++i)
  {
    sp_func_eval(mapping, &xis[i], &xjs[3*i]);
    xs[i].x = xjs[3*i];
    xs[i].y = xjs[3*i+1];
    xs[i].z = xjs[3*i+2];
  }
  // FIXME
}

void test_linear_block(void** state)
{
  block_t* block = block_new(1);
  assert_int_equal(1, block_order(block));
  assert_int_equal(8, block_num_points(block));
  for (int b = 0; b < 7; ++b)
    test_ball_mapping(block, b);
}

void test_quadratic_block(void** state)
{
  block_t* block = block_new(2);
  assert_int_equal(2, block_order(block));
  assert_int_equal(27, block_num_points(block));
  for (int b = 0; b < 7; ++b)
    test_ball_mapping(block, b);
}

void test_cubic_block(void** state)
{
  block_t* block = block_new(3);
  assert_int_equal(3, block_order(block));
  assert_int_equal(64, block_num_points(block));
  for (int b = 0; b < 7; ++b)
    test_ball_mapping(block, b);
}

void test_block_tags(void** state)
{
  block_t* block = block_new(1);
  block_add_tag(block, 0, "left");
  block_add_tag(block, 0, "left2");
  block_add_tag(block, 1, "right");
  int pos = 0;
  char* tag;

  bool result = block_next_tag(block, 0, &pos, &tag);
  assert_true(result);
  assert_true(!strcmp(tag, "left") || !strcmp(tag, "left2"));
  result = block_next_tag(block, 0, &pos, &tag);
  assert_true(result);
  assert_true(!strcmp(tag, "left2") || !strcmp(tag, "left"));
  result = block_next_tag(block, 0, &pos, &tag);
  assert_false(result);

  pos = 0;
  result = block_next_tag(block, 1, &pos, &tag);
  assert_true(result);
  assert_true(!strcmp(tag, "right"));
  result = block_next_tag(block, 1, &pos, &tag);
  assert_false(result);

  for (int f = 2; f < 6; ++f)
  {
    pos = 0;
    result = block_next_tag(block, 2, &pos, &tag);
    assert_false(result);
  }

  block_delete_tag(block, 0, "left2");
  pos = 0;
  result = block_next_tag(block, 0, &pos, &tag);
  assert_true(result);
  assert_true(!strcmp(tag, "left"));
  result = block_next_tag(block, 0, &pos, &tag);
  assert_false(result);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_linear_block),
    unit_test(test_quadratic_block),
    unit_test(test_cubic_block),
    unit_test(test_block_tags)
  };
  return run_tests(tests);
}
