#include "geometry/sphere.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
  point_t x;
  double r;
  normal_orient_t orient;
} sphere_t;

static void sphere_eval(void* ctx, point_t* x, double* result)
{
  sphere_t* s = ctx;
  double sign = (s->orient == OUTWARD_NORMAL) ? -1.0 : 1.0;
  result[0] = sign * (point_distance(x, &s->x) - s->r);
}

static void sphere_eval_gradient(void* ctx, point_t* x, double* result)
{
  sphere_t* s = ctx;
  double sign = (s->orient == OUTWARD_NORMAL) ? -1.0 : 1.0;
  double D = point_distance(x, &s->x);
  result[0] = sign * (x->x - s->x.x) / D;
  result[1] = sign * (x->x - s->x.x) / D;
  result[2] = sign * (x->x - s->x.x) / D;
}

sp_func_t* sphere_new(point_t* x, double r, normal_orient_t normal_orientation)
{
  // Set up a sphere signed distance function.
  sphere_t* s = malloc(sizeof(sphere_t));
  point_copy(&s->x, x);
  s->r = r;
  s->orient = normal_orientation;

  char sphere_str[1024];
  sprintf(sphere_str, "Sphere (x = (%g %g %g), r = %g)", 
          x->x, x->y, x->z, r);
  sp_vtable vtable = {.eval = sphere_eval, .dtor = free};
  sp_func_t* sphere = sp_func_new(sphere_str, s, vtable, SP_INHOMOGENEOUS, 1);

  // Register the gradient function.
  char sphere_grad_str[1024];
  sprintf(sphere_grad_str, "Sphere gradient (x = (%g %g %g), r = %g)", 
          x->x, x->y, x->z, r);
  sp_vtable vtable_g = {.eval = sphere_eval_gradient}; // Notice no dtor.
  sp_func_t* sphere_grad = sp_func_new(sphere_grad_str, s, vtable_g, SP_INHOMOGENEOUS, 3);
  sp_func_register_deriv(sphere, 1, sphere_grad);

  return sphere;
}

#ifdef __cplusplus
}
#endif

