// Copyright 2012-2013 Jeffrey Johnson.
// 
// This file is part of Polymec, and is licensed under the Apache License, 
// Version 2.0 (the "License"); you may not use this file except in 
// compliance with the License. You may may find the text of the license in 
// the LICENSE file at the top-level source directory, or obtain a copy of 
// it at
// 
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "geometry/create_boundary_generators.h"
#include "core/slist.h"
#include "core/unordered_map.h"
#include "core/kd_tree.h"

void create_boundary_generators(ptr_array_t* surface_points, 
                                ptr_array_t* surface_normals, 
                                ptr_array_t* surface_tags,
                                point_t** boundary_generators,
                                int* num_boundary_generators,
                                char*** tag_names,
                                int_array_t*** tags,
                                int* num_tags)
{
  ASSERT(surface_points->size >= 4); // surface must be closed!
  ASSERT(surface_points->size == surface_normals->size);
  ASSERT(surface_points->size == surface_tags->size);

  int num_surface_points = surface_points->size;

  // Compute the minimum distance from each surface point to its neighbors.
  double* h_min = malloc(sizeof(point_t) * num_surface_points);
  {
    // Dump the surface points into a kd-tree.
    point_t* surf_points = malloc(sizeof(point_t) * num_surface_points);
    for (int i = 0; i < num_surface_points; ++i)
      surf_points[i] = *((point_t*)surface_points->data[i]);
    kd_tree_t* tree = kd_tree_new(surf_points, num_surface_points);

    int neighbors[2];
    for (int i = 0; i < num_surface_points; ++i)
    {
      // Find the "nearest 2" points to the ith surface point--the first is 
      // the point itself, and the second is its nearest neighbor.
      kd_tree_nearest_n(tree, &surf_points[i], 2, neighbors);
printf("neighbors for point %d are %d and %d\n", i, neighbors[0], neighbors[1]);
      ASSERT(neighbors[0] == i);
      ASSERT(neighbors[1] >= 0);
      ASSERT(neighbors[1] < num_surface_points);
      h_min[i] = point_distance(&surf_points[i], &surf_points[neighbors[1]]);
    }

    // Clean up.
    kd_tree_free(tree);
    free(surf_points);
  }

  // Generate boundary points for each surface point based on how many 
  // surfaces it belongs to. 
  ptr_array_t* boundary_points = ptr_array_new();
  string_int_unordered_map_t* tag_indices = string_int_unordered_map_new();
  ptr_array_t* boundary_tags = ptr_array_new();
  for (int i = 0; i < num_surface_points; ++i)
  {
    // Add any tags from this point to the set of existing boundary tags.
    string_slist_t* tags = surface_tags->data[i];
    for (string_slist_node_t* t_iter = tags->front; t_iter != NULL; t_iter = t_iter->next)
      string_int_unordered_map_insert(tag_indices, t_iter->value, tag_indices->size);

    // Retrieve the surface point.
    point_t* x_surf = surface_points->data[i];

    // Retrieve the list of normal vectors for this surface point.
    ptr_slist_t* normal_list = surface_normals->data[i];
    int num_normals = normal_list->size;
    vector_t normals[num_normals];
    int n_offset = 0;
    for (ptr_slist_node_t* n_iter = normal_list->front; n_iter != NULL; n_iter = n_iter->next)
      normals[n_offset++] = *((vector_t*)n_iter->value);

    // For now, let's keep things relatively simple.
    if (num_normals > 2)
    {
      polymec_error("create_boundary_generators: Too many normal vectors (%d) for surface point %d at (%g, %g, %g)",
                    num_normals, i, x_surf->x, x_surf->y, x_surf->z);
    }

    // Create boundary points based on this list of normals.
    if (num_normals == 1)
    {
      // This point only belongs to one surface, so we create boundary points 
      // on either side of it.
      point_t* x_out = malloc(sizeof(point_t));
      x_out->x = x_surf->x + h_min[i]*normals[0].x;
      x_out->y = x_surf->y + h_min[i]*normals[0].y;
      x_out->z = x_surf->z + h_min[i]*normals[0].z;
      ptr_array_append_with_dtor(boundary_points, x_out, DTOR(free));

      point_t* x_in = malloc(sizeof(point_t));
      x_in->x = x_surf->x - h_min[i]*normals[0].x;
      x_in->y = x_surf->y - h_min[i]*normals[0].y;
      x_in->z = x_surf->z - h_min[i]*normals[0].z;
      ptr_array_append_with_dtor(boundary_points, x_in, DTOR(free));
    }
    else if (num_normals == 2)
    {
      // This point appears at the interface between two surfaces.
      // (Or so it seems.)
      ASSERT(vector_dot(&normals[0], &normals[1]) < 0.0);

      point_t* x1 = malloc(sizeof(point_t));
      x1->x = x_surf->x + h_min[i]*normals[0].x;
      x1->y = x_surf->y + h_min[i]*normals[0].y;
      x1->z = x_surf->z + h_min[i]*normals[0].z;
      ptr_array_append_with_dtor(boundary_points, x1, DTOR(free));

      point_t* x2 = malloc(sizeof(point_t));
      x2->x = x_surf->x - h_min[i]*normals[1].x;
      x2->y = x_surf->y - h_min[i]*normals[1].y;
      x2->z = x_surf->z - h_min[i]*normals[1].z;
      ptr_array_append_with_dtor(boundary_points, x2, DTOR(free));
    }

    // Tag the boundary point appropriately.
    ptr_array_append(boundary_tags, tags); // Borrowed ref to tags.
  }

  // Move the surface points into a contiguous array.
  *boundary_generators = malloc(sizeof(point_t) * boundary_points->size);
  *num_boundary_generators = boundary_points->size;
  for (int i = 0; i < boundary_points->size; ++i)
  {
    (*boundary_generators)[i] = *((point_t*)boundary_points->data[i]);
  }

  // Transcribe the tags.
  *tag_names = malloc(sizeof(char*) * tag_indices->size);
  *tags = malloc(sizeof(int_array_t**) * tag_indices->size);
  *num_tags = tag_indices->size;
  char* tag_name;
  int pos = 0, tag_index;
  while (string_int_unordered_map_next(tag_indices, &pos, &tag_name, &tag_index))
  {
    (*tag_names)[tag_index] = strdup(tag_name);
    (*tags)[tag_index] = int_array_new();
    for (int j = 0; j < *num_boundary_generators; ++j)
    {
      string_slist_t* tags_for_point = boundary_tags->data[j];
      for (string_slist_node_t* iter = tags_for_point->front; iter != NULL; iter = iter->next)
        int_array_append((*tags)[tag_index], tag_index);
    }
  }

  // Clean up.
  string_int_unordered_map_free(tag_indices);
  free(h_min);
}


