#include "core/constant_st_func.h"
#include "geometry/plane.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
  // Plane parameters.
  vector_t n;
  point_t x;

  // Basis vectors (for projections).
  vector_t e1, e2, e3;
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
  sp_vtable vtable = {.eval = plane_eval, .dtor = free};
  sp_func_t* plane = sp_func_new("Plane (uninitialized)", p, vtable, SP_INHOMOGENEOUS, 1);

  // Set up all the innards.
  plane_reset(plane, n, x);
  return plane;
}

void plane_reset(sp_func_t* plane, vector_t* n, point_t* x)
{
  plane_t* p = sp_func_context(plane);
  vector_copy(&p->n, n);
  vector_normalize(&p->n);
  point_copy(&p->x, x);

  char plane_str[1024];
  sprintf(plane_str, "Plane (n = (%g %g %g), x = (%g %g %g))", 
          n->x, n->y, n->z, x->x, x->y, x->z);
  sp_func_rename(plane, plane_str);

  // Register the negative of the normal as the derivative.
  double nn[3];
  nn[0] = -n->x, nn[1] = -n->y, nn[2] = -n->z; 
  sp_func_t* G = constant_sp_func_new(3, nn);
  sp_func_register_deriv(plane, 1, G);

  // Set up our basis vectors.
  vector_copy(&p->e3, &p->n);

  // Pick an arbitrary vector, e1, that is perpendicular to e3. One of these
  // should work.
  if (p->e3.x != 0.0)
  {
    p->e1.y = 1.0, p->e1.z = 1.0; 
    p->e1.x = -(p->e3.y + p->e3.z) / p->e3.x;
  }
  else if (p->e3.y != 0.0)
  {
    p->e1.x = 1.0, p->e1.z = 1.0; 
    p->e1.y = -(p->e3.x + p->e3.z) / p->e3.y;
  }
  else if (p->e3.z != 0.0)
  {
    p->e1.x = 1.0, p->e1.y = 1.0; 
    p->e1.z = -(p->e3.x + p->e3.y) / p->e3.z;
  }
  vector_normalize(&p->e1);

  // e2 = e3 x e1.
  vector_cross(&p->e3, &p->e1, &p->e2);
  ASSERT(vector_mag(&p->e2) > 1e-14);
}

void plane_project(sp_func_t* plane, point_t* x, double* eta, double* xi)
{
  plane_t* p = sp_func_context(plane);
  vector_t v = {.x = x->x - p->x.x, .y = x->y - p->x.y, .z = x->z - p->x.z};
  double voe3 = vector_dot(&v, &p->e3);
  vector_t v_perp = {.x = v.x - voe3 * p->e3.x, .y = v.y - voe3 * p->e3.y, .z = v.z - voe3 * p->e3.z};
  *eta = vector_dot(&v_perp, &p->e1);
  *xi = vector_dot(&v_perp, &p->e2);
}

double plane_intersect_with_line(sp_func_t* plane, point_t* x0, vector_t* t)
{
  plane_t* p = sp_func_context(plane);
  double not = vector_dot(&p->n, t);
  if (not == 0.0) // No intersection!
    return -FLT_MAX;
  else
  {
    double numer = p->n.x * (p->x.x - x0->x) +
                   p->n.y * (p->x.y - x0->y) +
                   p->n.z * (p->x.z - x0->z);
    return numer / not;
  }
}

#ifdef __cplusplus
}
#endif

