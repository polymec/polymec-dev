#include <gc/gc.h>
#include "geometry/cvt_gen_dist.h"
#include "geometry/scaled.h"

#ifdef __cplusplus
extern "C" {
#endif

struct cvt_gen_dist_t 
{
  char* name;
  void* context;
  double scale_factor;
  cvt_gen_dist_vtable vtable;
};

static void cvt_gen_dist_free(void* ctx, void* dummy)
{
  cvt_gen_dist_t* dist = (cvt_gen_dist_t*)ctx;
  free(dist->name);
  if ((dist->context != NULL) && (dist->vtable.dtor != NULL))
    dist->vtable.dtor(dist->context);
}

static void project_point_to_boundary(sp_func_t* boundary, point_t* x)
{
  ASSERT(sp_func_has_deriv(boundary, 1));
  double D, grad_D[3];
  sp_func_eval(boundary, x, &D);
  sp_func_eval_deriv(boundary, 1, x, grad_D);
  vector_t normal = {.x = -grad_D[0], .y = -grad_D[1], .z = -grad_D[2]};
  vector_normalize(&normal);
//printf("%g %g %g (%g, %g, %g, %g) ->", x->x, x->y, x->z, D, normal.x, normal.y, normal.z);
  x->x += D * normal.x;
  x->y += D * normal.y;
  x->z += D * normal.z;
//printf("%g %g %g\n", x->x, x->y, x->z);
}

static void reflect_point_across_boundary(sp_func_t* boundary, point_t* x)
{
  ASSERT(sp_func_has_deriv(boundary, 1));
  double D, grad_D[3];
  sp_func_eval(boundary, x, &D);
  sp_func_eval_deriv(boundary, 1, x, grad_D);
  vector_t normal = {.x = -grad_D[0], .y = -grad_D[1], .z = -grad_D[2]};
  vector_normalize(&normal);
  x->x += 2.0 * D * normal.x;
  x->y += 2.0 * D * normal.y;
  x->z += 2.0 * D * normal.z;
}

cvt_gen_dist_t* cvt_gen_dist_new(const char* name, void* context, cvt_gen_dist_vtable vtable)
{
  ASSERT(vtable.iterate != NULL);

  cvt_gen_dist_t* dist = GC_MALLOC(sizeof(cvt_gen_dist_t));
  dist->name = strdup(name);
  dist->context = context;
  dist->vtable = vtable;
  dist->scale_factor = 1.0;
  GC_register_finalizer(dist, cvt_gen_dist_free, dist, NULL, NULL);
  return dist;
}

const char* cvt_gen_dist_name(cvt_gen_dist_t* dist)
{
  return (const char*)dist->name;
}

void cvt_gen_dist_set_safety_buffer(cvt_gen_dist_t* dist, double factor)
{
  ASSERT(factor >= 0.0);
  ASSERT(factor < 1.0);
  dist->scale_factor = 1.0 - factor;
}

void cvt_gen_dist_iterate(cvt_gen_dist_t* dist, 
                          sp_func_t* density,
                          sp_func_t* boundary,
                          bbox_t* bounding_box,
                          point_t* interior_points, 
                          int num_interior_points)
{
  cvt_gen_dist_iterate_with_boundary_points(dist, density, boundary, 
    bounding_box, interior_points, num_interior_points, NULL, 0);
}

void cvt_gen_dist_iterate_with_boundary_points(cvt_gen_dist_t* dist, 
                                               sp_func_t* density,
                                               sp_func_t* boundary,
                                               bbox_t* bounding_box,
                                               point_t* interior_points, 
                                               int num_interior_points,
                                               point_t* boundary_points,
                                               int num_boundary_points)
{
  sp_func_t* interior = boundary;
  if ((dist->scale_factor < 1.0) && (boundary != NULL))
    interior = scaled_new(boundary, dist->scale_factor);

  // If necessary, move interior points into the interior region.
  if ((num_interior_points > 0) && (interior != NULL))
  {
    for (int i = 0; i < num_interior_points; ++i)
    {
      double D;
      sp_func_eval(interior, &interior_points[i], &D);
      if (D >= 0.0)
      {
        // We attempt to put the point inside the domain by reflecting it 
        // across the boundary. If the boundary is some composite like an 
        // intersection, we may have to do more than one reflection.
        static const int max_reflections = 100;
        for (int j = 0; j < max_reflections; ++j) 
        {
          reflect_point_across_boundary(interior, &interior_points[i]);
          sp_func_eval(interior, &interior_points[i], &D);
          if (D < 0.0) break;
        }
        if (D >= 0.0)
          polymec_error("cvt_gen_dist_iterate: reflection of outside point yielded\n another outside point. Boundary is not a signed distance function!");
      }
    }
  }

  // Project the boundary points to the boundary.
  if ((num_boundary_points > 0) && (boundary != NULL))
  {
    for (int i = 0; i < num_boundary_points; ++i)
    {
      project_point_to_boundary(boundary, &boundary_points[i]);
      double D;
      sp_func_eval(boundary, &boundary_points[i], &D);
      if (fabs(D) > 1e-12)
        polymec_error("cvt_gen_dist_iterate: boundary projection yielded a non-zero distance (%g).\n Boundary is not a signed distance function!", D);
    }
  }

  dist->vtable.iterate(dist->context, density, bounding_box,
                       interior_points, num_interior_points, interior,
                       boundary_points, num_boundary_points, boundary);
}

#ifdef __cplusplus
}
#endif

