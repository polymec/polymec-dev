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

#if 0
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
#endif

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

static inline void swap_point(point_t* points, int i, int j)
{
  point_t temp;
  point_copy(&temp, &points[i]);
  point_copy(&points[i], &points[j]);
  point_copy(&points[j], &temp);
}

void cvt_gen_dist_iterate(cvt_gen_dist_t* dist, 
                          sp_func_t* density,
                          sp_func_t* boundary,
                          bbox_t* bounding_box,
                          point_t* points, 
                          int num_points,
                          int* num_boundary_points)
{
  ASSERT(num_points > 0);
  ASSERT(points != NULL);
  ASSERT(num_boundary_points != NULL);

  // If necessary, project points back to the boundary.
  if (boundary != NULL)
  {
    for (int i = 0; i < num_points; ++i)
    {
      double D;
      sp_func_eval(boundary, &points[i], &D);
      if (D >= 0.0)
      {
#if 0
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
#endif
        project_point_to_boundary(boundary, &points[i]);
        double D;
        sp_func_eval(boundary, &points[i], &D);
        if (fabs(D) > 1e-12)
          polymec_error("cvt_gen_dist_iterate: boundary projection yielded a non-zero distance (%g).\n Boundary is not a signed distance function!", D);
      }
    }
  }

  // Perform the iteration.
  dist->vtable.iterate(dist->context, density, boundary, bounding_box, 
                       points, num_points);

  // Check that all points are within the domain or on the boundary, and 
  // move the boundary points to the end of the list.
  *num_boundary_points = 0;
  if (boundary != NULL)
  {
    for (int i = 0; i < num_points; ++i)
    {
      double D;
      sp_func_eval(boundary, &points[i], &D);
      if (D > 0.0)
      {
        polymec_error("cvt_gen_dist_iterate: Point %d lies outside the domain (D = %g).", i, D);
      }
      else if (fabs(D) < 1e-12)
      {
        // This point is on the boundary, so we move it to the end of the list.
        ++(*num_boundary_points);
        swap_point(points, i, num_points - *num_boundary_points);
      }
    }
  }
}

#ifdef __cplusplus
}
#endif

