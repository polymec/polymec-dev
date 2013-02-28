#include "geometry/ball_mapping.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
  point_t x0;
  double r, l;
  int index;
} ball_mapping_t;

static void ball_map0(void* context, point_t* xi, double* x)
{
}

static void ball_map1(void* context, point_t* xi, double* x)
{
}

static void ball_map2(void* context, point_t* xi, double* x)
{
}

static void ball_map3(void* context, point_t* xi, double* x)
{
}

static void ball_map4(void* context, point_t* xi, double* x)
{
}

static void ball_map5(void* context, point_t* xi, double* x)
{
}

static void ball_map6(void* context, point_t* xi, double* x)
{
}

static void ball_dtor(void* context)
{
  ball_mapping_t* ball = context;
  free(ball);
}

sp_func_t* ball_mapping_new(point_t* x0, double r, double l, int block_index)
{
  ASSERT(block_index >= 0);
  ASSERT(block_index <= 6);
  sp_vtable vtable = {.dtor = ball_dtor};
  switch (block_index)
  {
    case 0:
      vtable.eval = ball_map0;
      break;
    case 1:
      vtable.eval = ball_map1;
      break;
    case 2:
      vtable.eval = ball_map2;
      break;
    case 3:
      vtable.eval = ball_map3;
      break;
    case 4:
      vtable.eval = ball_map4;
      break;
    case 5:
      vtable.eval = ball_map5;
      break;
    case 6:
      vtable.eval = ball_map6;
      break;
  }
  char name[1024];
  snprintf(name, 1024, "Multi-block ball mapping (block %d)", block_index);
  ball_mapping_t* m = malloc(sizeof(ball_mapping_t));
  m->r = r;
  m->x0.x = x0->x;
  m->x0.y = x0->y;
  m->x0.z = x0->z;
  m->index = block_index;
  return sp_func_new(name, (void*)m, vtable, SP_INHOMOGENEOUS, 3);
}

#ifdef __cplusplus
}
#endif

