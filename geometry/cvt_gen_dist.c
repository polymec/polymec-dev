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
  dist->vtable.iterate(dist->context, density, bounding_box,
                       interior_points, num_interior_points, interior,
                       boundary_points, num_boundary_points, boundary);
}

#ifdef __cplusplus
}
#endif

