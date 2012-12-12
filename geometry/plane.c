#include "core/constant_st_func.h"
#include "geometry/plane.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
  vector_t n;
  point_t x;
} plane_t;

static void plane_eval(void* ctx, point_t* x, double* result)
{
  plane_t* p = ctx;
  vector_t D;
  point_displacement(&p->x, x, &D);
  result[0] = vector_dot(&p->n, &D);
}

sp_func_t* plane_new(vector_t* n, point_t* x)
{
  // Set up a plane signed distance function.
  plane_t* p = malloc(sizeof(plane_t));
  vector_copy(&p->n, n);
  point_copy(&p->x, x);

  char plane_str[1024];
  sprintf(plane_str, "Plane (n = (%g %g %g), x = (%g %g %g))", 
          n->x, n->y, n->z, x->x, x->y, x->z);
  sp_vtable vtable = {.eval = plane_eval, .dtor = free};
  sp_func_t* plane = sp_func_new(plane_str, p, vtable, SP_INHOMOGENEOUS, 1);

  // Register the negative of the normal as the derivative.
  double nn[3];
  nn[0] = -n->x, nn[1] = -n->y, nn[2] = -n->z; 
  sp_func_t* G = constant_sp_func_new(3, nn);
  sp_func_register_deriv(plane, 1, G);

  return plane;
}

#ifdef __cplusplus
}
#endif

