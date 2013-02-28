#include "geometry/cylinder_mapping.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
  point_t x0;
  vector_t axis;
  double r, l;
  int index;
} cylinder_mapping_t;

static void cylinder_map0(void* context, point_t* xi, double* x)
{
}

static void cylinder_map1(void* context, point_t* xi, double* x)
{
}

static void cylinder_map2(void* context, point_t* xi, double* x)
{
}

static void cylinder_map3(void* context, point_t* xi, double* x)
{
}

static void cylinder_map4(void* context, point_t* xi, double* x)
{
}

static void cylinder_dtor(void* context)
{
  cylinder_mapping_t* ball = context;
  free(ball);
}

sp_func_t* cylinder_mapping_new(point_t* x0, vector_t* axis, double r, double l, int block_index)
{
  ASSERT(block_index >= 0);
  ASSERT(block_index <= 4);
  sp_vtable vtable = {.dtor = cylinder_dtor};
  switch (block_index)
  {
    case 0:
      vtable.eval = cylinder_map0;
      break;
    case 1:
      vtable.eval = cylinder_map1;
      break;
    case 2:
      vtable.eval = cylinder_map2;
      break;
    case 3:
      vtable.eval = cylinder_map3;
      break;
    case 4:
      vtable.eval = cylinder_map4;
      break;
  }
  char name[1024];
  snprintf(name, 1024, "Multi-block cylinder mapping (block %d)", block_index);
  cylinder_mapping_t* m = malloc(sizeof(cylinder_mapping_t));
  m->r = r;
  m->l = l;
  m->x0.x = x0->x;
  m->x0.y = x0->y;
  m->x0.z = x0->z;
  m->axis.x = axis->x;
  m->axis.y = axis->y;
  m->axis.z = axis->z;
  m->index = block_index;
  return sp_func_new(name, (void*)m, vtable, SP_INHOMOGENEOUS, 3);
}

#ifdef __cplusplus
}
#endif

