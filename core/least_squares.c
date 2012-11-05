#include <stdlib.h>
#include <gc/gc.h>
#include "core/least_squares.h"

#ifdef __cplusplus
extern "C" {
#endif

struct multi_index_t 
{
  int p, x_order, y_order, z_order;
  int offset; // Offset in flattened index space.
};

static void multi_index_free(void* ctx, void* dummy)
{
  multi_index_t* m = (multi_index_t*)ctx;
  free(m);
}

multi_index_t* multi_index_new(int p)
{
  ASSERT(p >= 0);
  ASSERT(p < 4);
  multi_index_t* m = GC_MALLOC(sizeof(multi_index_t));
  m->p = p;
  m->x_order = m->y_order = m->z_order = 0;
  m->offset = 0;
  GC_register_finalizer(m, &multi_index_free, m, NULL, NULL);
  return m;
}

static int multi_index_sizes[] = {1, 4, 10, 18};

bool multi_index_next(multi_index_t* m, int* x_order, int* y_order, int* z_order)
{
  ASSERT(m->p >= 0);
  ASSERT(m->p < 4);
  if (m->p == 0)
  {
    *x_order = m->x_order;
    *y_order = m->y_order;
    *z_order = m->z_order;
  }
  if (m->p == 1)
  {
    static const multi_index_t multi_index[] = 
    {
      {.p = 1, .x_order = 0, .y_order = 0, .z_order = 0, .offset = 0}, 
      {.p = 1, .x_order = 1, .y_order = 0, .z_order = 0, .offset = 1}, 
      {.p = 1, .x_order = 0, .y_order = 1, .z_order = 0, .offset = 2}, 
      {.p = 1, .x_order = 0, .y_order = 0, .z_order = 1, .offset = 3},
      {.p = -1}
    };
    *x_order = m->x_order;
    *y_order = m->y_order;
    *z_order = m->z_order;
    *m = multi_index[m->offset+1];
  }
  else if (m->p == 2)
  {
    static const multi_index_t multi_index[] = 
    {
      {.p = 2, .x_order = 0, .y_order = 0, .z_order = 0, .offset = 0}, 
      {.p = 2, .x_order = 1, .y_order = 0, .z_order = 0, .offset = 1}, 
      {.p = 2, .x_order = 0, .y_order = 1, .z_order = 0, .offset = 2}, 
      {.p = 2, .x_order = 0, .y_order = 0, .z_order = 1, .offset = 3},
      {.p = 2, .x_order = 2, .y_order = 0, .z_order = 0, .offset = 4},
      {.p = 2, .x_order = 1, .y_order = 1, .z_order = 0, .offset = 5},
      {.p = 2, .x_order = 1, .y_order = 0, .z_order = 1, .offset = 6},
      {.p = 2, .x_order = 0, .y_order = 2, .z_order = 0, .offset = 7},
      {.p = 2, .x_order = 0, .y_order = 1, .z_order = 1, .offset = 8},
      {.p = 2, .x_order = 0, .y_order = 0, .z_order = 2, .offset = 9},
      {.p = -1}
    };
    *x_order = m->x_order;
    *y_order = m->y_order;
    *z_order = m->z_order;
    *m = multi_index[m->offset+1];
  }
  else if (m->p == 3)
  {
    static const multi_index_t multi_index[] = 
    {
      {.p = 3, .x_order = 0, .y_order = 0, .z_order = 0, .offset = 0}, 
      {.p = 3, .x_order = 1, .y_order = 0, .z_order = 0, .offset = 1}, 
      {.p = 3, .x_order = 0, .y_order = 1, .z_order = 0, .offset = 2}, 
      {.p = 3, .x_order = 0, .y_order = 0, .z_order = 1, .offset = 3},
      {.p = 3, .x_order = 2, .y_order = 0, .z_order = 0, .offset = 4},
      {.p = 3, .x_order = 1, .y_order = 1, .z_order = 0, .offset = 5},
      {.p = 3, .x_order = 1, .y_order = 0, .z_order = 1, .offset = 6},
      {.p = 3, .x_order = 0, .y_order = 2, .z_order = 0, .offset = 7},
      {.p = 3, .x_order = 0, .y_order = 1, .z_order = 1, .offset = 8},
      {.p = 3, .x_order = 0, .y_order = 0, .z_order = 2, .offset = 9},
      {.p = 3, .x_order = 3, .y_order = 0, .z_order = 0, .offset = 10},
      {.p = 3, .x_order = 2, .y_order = 1, .z_order = 0, .offset = 11},
      {.p = 3, .x_order = 2, .y_order = 0, .z_order = 1, .offset = 12},
      {.p = 3, .x_order = 1, .y_order = 2, .z_order = 0, .offset = 13},
      {.p = 3, .x_order = 1, .y_order = 1, .z_order = 1, .offset = 14},
      {.p = 3, .x_order = 0, .y_order = 3, .z_order = 0, .offset = 15},
      {.p = 3, .x_order = 0, .y_order = 2, .z_order = 1, .offset = 16},
      {.p = 3, .x_order = 0, .y_order = 1, .z_order = 2, .offset = 17},
      {.p = -1}
    };
    *x_order = m->x_order;
    *y_order = m->y_order;
    *z_order = m->z_order;
    *m = multi_index[m->offset+1];
  }
  return ((m->offset+1) < multi_index_sizes[m->p]);
}

void multi_index_reset(multi_index_t* m)
{
  m->offset = 0;
  m->x_order = m->y_order = m->z_order = 0;
}

int multi_index_order(multi_index_t* m)
{
  return m->p;
}

int multi_index_size(multi_index_t* m)
{
  return multi_index_sizes[m->p];
}

void compute_moment_matrix(int p, point_t* points, int num_points, double* A)
{
  multi_index_t* m = multi_index_new(p);
  int x, y, z, offset = 0;
  int size = multi_index_size(m);
  for (int i = 0; i < num_points; ++i)
  {
    while (multi_index_next(m, &x, &y, &z))
      A[i*size + offset] = pow(points[i].x, x)*pow(points[i].y, y)*pow(points[i].z, z);
  }
}

#ifdef __cplusplus
}
#endif

