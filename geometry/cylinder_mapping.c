#include "geometry/cylinder_mapping.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
  point_t x0;
  double r, l;
  int index;
} cylinder_mapping_t;

static void cylinder_map0(void* context, point_t* xi, double* x)
{
  // Inner cube.
  cylinder_mapping_t* cyl = context;
  x[0] = cyl->x0.x + 0.5*cyl->r * (xi->x - 0.5);
  x[1] = cyl->x0.y + 0.5*cyl->r * (xi->y - 0.5);
  x[2] = cyl->x0.z + cyl->l * (xi->z - 0.5);
}

static void cylinder_map1(void* context, point_t* xi, double* x)
{
  // -x cube.
  cylinder_mapping_t* cyl = context;

  double Xi = M_PI/2 * (0.5 - xi->y);
  double phi = -M_PI + Xi;
  double r0 = 0.5*cyl->r / cos(Xi); // Closest approach
  double r = cyl->r - xi->x * (cyl->r - r0);

  x[0] = cyl->x0.x + r * cos(phi);
  x[1] = cyl->x0.y + r * sin(phi);
  x[2] = cyl->x0.z + cyl->l * (xi->z - 0.5);
}

static void cylinder_map2(void* context, point_t* xi, double* x)
{
  // +x cube.
  cylinder_mapping_t* cyl = context;

  double Xi = M_PI/2 * (xi->y - 0.5);
  double phi = Xi;
  double r0 = 0.5*cyl->r / cos(Xi); // Closest approach
  double r = xi->x * (cyl->r - r0) - cyl->r;

  x[0] = cyl->x0.x + r * cos(phi);
  x[1] = cyl->x0.y + r * sin(phi);
  x[2] = cyl->x0.z + cyl->l * (xi->z - 0.5);
}

static void cylinder_map3(void* context, point_t* xi, double* x)
{
  // -y cube.
  cylinder_mapping_t* cyl = context;

  double Xi = M_PI/2 * (xi->x - 0.5);
  double phi = -M_PI/2 + Xi;
  double r0 = 0.5*cyl->r / cos(Xi); // Closest approach
  double r = cyl->r - xi->y * (cyl->r - r0);

  x[0] = cyl->x0.x + r * cos(phi);
  x[1] = cyl->x0.y + r * sin(phi);
  x[2] = cyl->x0.z + cyl->l * (xi->z - 0.5);
}

static void cylinder_map4(void* context, point_t* xi, double* x)
{
  // +y cube.
  cylinder_mapping_t* cyl = context;

  double Xi = M_PI/2 * (0.5 - xi->x);
  double phi = M_PI/2 + Xi;
  double r0 = 0.5*cyl->r / cos(Xi); // Closest approach
  double r = xi->y * (cyl->r - r0) - cyl->r;

  x[0] = cyl->x0.x + r * cos(phi);
  x[1] = cyl->x0.y + r * sin(phi);
  x[2] = cyl->x0.z + cyl->l * (xi->z - 0.5);
}

static void cylinder_dtor(void* context)
{
  cylinder_mapping_t* ball = context;
  free(ball);
}

sp_func_t* cylinder_mapping_new(point_t* x0, double r, double l, int block_index)
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
  m->index = block_index;
  return sp_func_new(name, (void*)m, vtable, SP_INHOMOGENEOUS, 3);
}

#ifdef __cplusplus
}
#endif

