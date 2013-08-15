#include "geometry/create_boundary_generators.h"

void create_boundary_generators(ptr_array_t* surface_points, 
                                ptr_array_t* surface_normals, 
                                ptr_array_t* surface_tags,
                                point_t** boundary_generators,
                                int* num_boundary_generators,
                                char*** tag_names,
                                int_array_t*** tags,
                                int* num_tags)
{
  ASSERT(surface_points->size > 0);
  ASSERT(surface_points->size == surface_normals->size);
  ASSERT(surface_points->size == surface_tags->size);

}


