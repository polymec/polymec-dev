#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmockery.h"
#include "geometry/ball_mapping.h"

void plot_points(point_t* points, int num_points, const char* filename)
{
  FILE* fd = fopen(filename, "w");
  fprintf(fd, "# x y z\n");
  for (int i = 0; i < num_points; ++i)
    fprintf(fd, "%g %g %g\n", points[i].x, points[i].y, points[i].z);
  fclose(fd);
}

static const double L = 1.0;
static const int num_test_points = 27;
static point_t test_xs[7][27] = {
  // block 0
  {{.x = -0.5, .y = - 0.5, .z = -0.5},    // xi = (0, 0, 0)
   {.x = -0.5, .y = - 0.5, .z = 0.},      // xi = (0, 0, 0.5)
   {.x = -0.5, .y = - 0.5, .z = 0.5},     // xi = (0, 0, 1)
   {.x = -0.5, .y = 0., .z = -0.5},       // xi = (0, 0.5, 0)
   {.x = -0.5, .y = 0., .z = 0.},         // xi = (0, 0.5, 0.5)
   {.x = -0.5, .y = 0., .z = 0.5},        // xi = (0, 0.5, 1)
   {.x = -0.5, .y = 0.5, .z = -0.5},      // xi = (0, 1, 0)
   {.x = -0.5, .y = 0.5, .z = 0.},        // xi = (0, 1, 0.5)
   {.x = -0.5, .y = 0.5, .z = 0.5},       // xi = (0, 1, 1)
   {.x = 0., .y = -0.5, .z = -0.5},       // xi = (0.5, 0, 0)
   {.x = 0., .y = -0.5, .z = 0.},         // xi = (0.5, 0, 0.5)
   {.x = 0., .y = -0.5, .z = 0.5},        // xi = (0.5, 0, 1)
   {.x = 0., .y = 0., .z = -0.5},         // xi = (0.5, 0.5, 0)
   {.x = 0., .y = 0., .z = 0.},           // xi = (0.5, 0.5, 0.5)
   {.x = 0., .y = 0., .z = 0.5},          // xi = (0.5, 0.5, 1)
   {.x = 0., .y = 0.5, .z = -0.5},        // xi = (0.5, 1, 0)
   {.x = 0., .y = 0.5, .z = 0.},          // xi = (0.5, 1, 0.5)
   {.x = 0., .y = 0.5, .z = 0.5},         // xi = (0.5, 1, 1)
   {.x = 0.5, .y = -0.5, .z = -0.5},      // xi = (1, 0, 0)
   {.x = 0.5, .y = -0.5, .z = 0.},        // xi = (1, 0, 0.5)
   {.x = 0.5, .y = -0.5, .z = 0.5},       // xi = (1, 0, 1)
   {.x = 0.5, .y = 0., .z = -0.5},        // xi = (1, 0.5, 0)
   {.x = 0.5, .y = 0., .z = 0.},          // xi = (1, 0.5, 0.5)
   {.x = 0.5, .y = 0., .z = 0.5},         // xi = (1, 0.5, 1)
   {.x = 0.5, .y = 0.5, .z = -0.5},       // xi = (1, 1, 0)
   {.x = 0.5, .y = 0.5, .z = 0.},         // xi = (1, 1, 0.5)
   {.x = 0.5, .y = 0.5, .z = 0.5}},       // xi = (1, 1, 1)
  // block 1
};

static void test_block(int block_index)
{
  point_t zero = {.x = 0.0, .y = 0.0, .z = 0.0};
  sp_func_t* mapping = ball_mapping_new(&zero, 1.0, 0.5, block_index);
  point_t test_points[num_test_points], xis[num_test_points];

  // Construct the test points in logical space.
  int offset = 0;
  for (int i = 0; i < 3; ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      for (int k = 0; k < 3; ++k, ++offset)
      {
        xis[offset].x = 0.5*i;
        xis[offset].y = 0.5*j;
        xis[offset].z = 0.5*k;
      }
    }
  }
  ASSERT(offset == num_test_points);

  for (int i = 0; i < num_test_points; ++i)
  {
    double X[3];
    sp_func_eval(mapping, &xis[i], X);
    test_points[i].x = X[0];
    test_points[i].y = X[1];
    test_points[i].z = X[2];

//    double error = point_distance(&test_points[i], &test_xs[block_index][i]);
//    assert_true(error < 1e-12);
  }

  char fn[1024];
  snprintf(fn, 1024, "points.%d", block_index);
  plot_points(test_points, 27, fn);
}

void test_block0(void** state)
{
  test_block(0);
}

void test_block1(void** state)
{
  test_block(1);
}

void test_block2(void** state)
{
  test_block(2);
}

void test_block3(void** state)
{
  test_block(3);
}

void test_block4(void** state)
{
  test_block(4);
}

void test_block5(void** state)
{
  test_block(5);
}

void test_block6(void** state)
{
  test_block(6);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_block0),
    unit_test(test_block1),
    unit_test(test_block2),
    unit_test(test_block3),
    unit_test(test_block4),
    unit_test(test_block5),
    unit_test(test_block6)
  };
  return run_tests(tests);
}
