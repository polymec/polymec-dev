#include "geometry/create_cvt.h"
#include "geometry/create_voronoi_mesh.h"

struct cvt_iterator_t 
{
  char* name;
  void* context;
  cvt_iterator_vtable vtable;
};

cvt_iterator_t* cvt_iterator_new(const char* name, void* context, cvt_iterator_vtable vtable)
{
  ASSERT(vtable.move_points != NULL);
  ASSERT(vtable.is_finished != NULL);
  cvt_iterator_t* iter = malloc(sizeof(cvt_iterator_t));
  iter->name = strdup(name);
  iter->context = context;
  iter->vtable = vtable;
  return iter;
}

static void cvt_iterator_free(cvt_iterator_t* cvt_iter)
{
  if ((cvt_iter->context != NULL) && (cvt_iter->vtable.dtor != NULL))
    (cvt_iter->vtable.dtor)(cvt_iter->context);
  free(cvt_iter->name);
  free(cvt_iter);
}

mesh_t* create_cvt(point_t* stationary_generators, int num_stationary_generators, 
                   point_t* mobile_generators, int num_mobile_generators,
                   cvt_iterator_t* cvt_iter)
{
  ASSERT(num_stationary_generators >= 0);
  ASSERT(num_mobile_generators > 0);
  ASSERT(cvt_iter != NULL);

  // Initialize our iterator if needed.
  if (cvt_iter->vtable.init != NULL)
  {
    cvt_iter->vtable.init(cvt_iter->context, stationary_generators, num_stationary_generators,
                          mobile_generators, num_mobile_generators);
  }

  // Create an initial tessellation from all the points.
  int num_generators = num_stationary_generators + num_mobile_generators;
  point_t* all_generators = malloc(sizeof(point_t) * num_generators);
  point_t* my_mobile_generators = malloc(sizeof(point_t) * num_mobile_generators);
  memcpy(all_generators, stationary_generators, sizeof(point_t) * num_stationary_generators);
  memcpy(&all_generators[num_stationary_generators], mobile_generators, sizeof(point_t) * num_mobile_generators);
  memcpy(my_mobile_generators, mobile_generators, sizeof(point_t) * num_mobile_generators);
  mesh_t* mesh = create_voronoi_mesh(all_generators, num_generators, NULL, 0);

  // Iterate till we're done.
  int iteration = 0;
  while (!cvt_iter->vtable.is_finished(cvt_iter->context, mesh, iteration))
  {
    // Move the mobile points.
    cvt_iter->vtable.move_points(cvt_iter->context, my_mobile_generators, num_mobile_generators);

    // Destroy and recreate the mesh.
    mesh_free(mesh);
    memcpy(all_generators, stationary_generators, sizeof(point_t) * num_stationary_generators);
    memcpy(&all_generators[num_stationary_generators], mobile_generators, sizeof(point_t) * num_mobile_generators);
    mesh = create_voronoi_mesh(all_generators, num_generators, NULL, 0);
    ++iteration;
  }

  // Clean up.
  free(my_mobile_generators);
  free(all_generators);
  cvt_iterator_free(cvt_iter);

  return mesh;
}


