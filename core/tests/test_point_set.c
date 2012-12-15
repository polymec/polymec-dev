#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmockery.h"
#include "core/point_set.h"

void test_construct(void** state) 
{ 
  // Create a point set containing 100 random points.
  int N = 100;
  bbox_t bounding_box = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  point_set_t* pset = point_set_new();
  for (int i = 0; i < N; ++i)
  {
    point_t p;
    point_randomize(&p, random, &bounding_box);
    point_set_insert(pset, &p, i);
  }

  assert_int_equal(100, point_set_size(pset));
  point_set_clear(pset);
  assert_int_equal(0, point_set_size(pset));

  point_set_free(pset);
}

void test_find_nearest(void** state) 
{ 
  // Create a point set containing 100 random points.
  int N = 100;
  bbox_t bounding_box = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
  point_set_t* pset = point_set_new();
  point_t points[N];
  for (int i = 0; i < N; ++i)
  {
    point_randomize(&points[i], random, &bounding_box);
    point_set_insert(pset, &points[i], i);
  }

  // Now do some nearest point queries and check the answers.
  for (int i = 0; i < 10; ++i) // 10 queries.
  {
    // Pick a random point p.
    point_t p;
    point_randomize(&p, random, &bounding_box);

    // Which point in the set is p closest to?
    int j = point_set_nearest(pset, &p); // Point set says: "j."

    // Do a linear search for the closest point.
    double min_dist = FLT_MAX;
    int kk = -1;
    for (int k = 0; k < N; ++k)
    {
      double dist = point_distance(&p, &points[k]);
      if (dist < min_dist)
      {
        kk = k;
        min_dist = dist;
      }
    }

    // The two queries must agree with each other.
    assert_true(j == kk);
  }

  point_set_free(pset);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_construct),
    unit_test(test_find_nearest)
  };
  return run_tests(tests);
}
