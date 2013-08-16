#include "geometry/create_boundary_generators.h"
#include "core/unordered_map.h"

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

  // Generate boundary points for each surface point based on how many 
  // surfaces it belongs to. 
  int num_surface_points = surface_points->size;
  ptr_array_t* boundary_points = ptr_array_new();
  string_array_t* boundary_tag_names = string_array_new();
  string_int_unordered_map_t* tag_indices = string_int_unordered_map_new();
  for (int i = 0; i < num_surface_points; ++i)
  {
  }

  // Move the surface points into a contiguous array.
  *boundary_generators = malloc(sizeof(point_t) * boundary_points->size);
  *num_boundary_generators = boundary_points->size;
  for (int i = 0; i < boundary_points->size; ++i)
  {
    (*boundary_generators)[i] = *((point_t*)boundary_points->data[i]);
  }

  // Transcribe the tags.
  *tag_names = malloc(sizeof(char*) * boundary_tag_names->size);
  *tags = malloc(sizeof(int_array_t**) * boundary_tag_names->size);
  *num_tags = boundary_tag_names->size;
  for (int i = 0; i < boundary_tag_names->size; ++i)
  {
    (*tag_names)[i] = strdup(boundary_tag_names->data[i]);
    (*tags)[i] = int_array_new();
    for (int j = 0; j < num_surface_points; ++j)
    {
      string_array_t* tags_for_point = surface_tags->data[j];
      for (int k = 0; k < tags_for_point->size; ++k)
      {
        int tag_index = *string_int_unordered_map_get(tag_indices, tags_for_point->data[k]);
        int_array_append((*tags)[i], tag_index);
      }
    }
  }

  // Clean up.
  string_int_unordered_map_free(tag_indices);
  string_array_free(boundary_tag_names);
}


