#include "geometry/cylinder.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
  vector_t d;
  point_t x;
  double r;
  normal_orient_t orient;
} cyl_t;

static void cyl_eval(void* ctx, point_t* x, double* result)
{
  cyl_t* c = ctx;
  double sign = (c->orient == OUTWARD_NORMAL) ? -1.0 : 1.0;
  vector_t disp, rad_disp;
  point_displacement(&c->x, x, &disp);
  vector_cross(&disp, &c->d, &rad_disp); // Right??
  result[0] = sign * (vector_mag(&rad_disp) - c->r);
}

static void cyl_eval_gradient(void* ctx, point_t* x, double* result)
{
  cyl_t* c = ctx;
  double sign = (c->orient == OUTWARD_NORMAL) ? -1.0 : 1.0;
  vector_t disp, rad_disp;
  point_displacement(&c->x, x, &disp);
  vector_cross(&disp, &c->d, &rad_disp); // Right??
  double D = vector_mag(&rad_disp);
  if (D == 0.0)
  {
    result[0] = result[1] = result[2] = 0.0;
  }
  else
  {
    result[0] = sign * rad_disp.x / D;
    result[1] = sign * rad_disp.y / D;
    result[2] = sign * rad_disp.z / D;
  }
}

sp_func_t* cylinder_new(vector_t* d, point_t* x, double r, normal_orient_t normal_orientation)
{
  // Set up a cylinder signed distance function.
  cyl_t* c = malloc(sizeof(cyl_t));
  vector_copy(&c->d, d);
  vector_normalize(&c->d);
  point_copy(&c->x, x);
  c->r = r;
  c->orient = normal_orientation;

  char cyl_str[1024];
  sprintf(cyl_str, "Cylinder (d = (%g, %g, %g), x = (%g %g %g), r = %g)", 
          d->x, d->y, d->z, x->x, x->y, x->z, r);
  sp_vtable vtable = {.eval = cyl_eval, .dtor = free};
  sp_func_t* cyl = sp_func_new(cyl_str, c, vtable, SP_INHOMOGENEOUS, 1);

  // Register the gradient function.
  char cyl_grad_str[1024];
  sprintf(cyl_grad_str, "Cylinder gradient (d = (%g, %g, %g), x = (%g %g %g), r = %g)", 
          d->x, d->y, d->z, x->x, x->y, x->z, r);
  sp_vtable vtable_g = {.eval = cyl_eval_gradient}; // Notice no dtor.
  sp_func_t* cyl_grad = sp_func_new(cyl_grad_str, c, vtable_g, SP_INHOMOGENEOUS, 3);
  sp_func_register_deriv(cyl, 1, cyl_grad);

  return cyl;
}

#ifdef __cplusplus
}
#endif

