// Copyright (c) 2012-2013, Jeffrey N. Johnson
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this 
// list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice, 
// this list of conditions and the following disclaimer in the documentation 
// and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


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
  double* h_min = malloc(sizeof(double) * num_surface_points);
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
      // FIXME: Serious memory error within here. We work around it for 
      // FIXME: the moment by freshing the tree pointer, which is 
      // FIXME: corrupted. ICK!
      kd_tree_t* tree_p = tree;
      kd_tree_nearest_n(tree, &surf_points[i], 2, neighbors);
      tree = tree_p;
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
    (*tag_names)[tag_index] = string_dup(tag_name);
    (*tags)[tag_index] = int_array_new();
    for (int j = 0; j < *num_boundary_generators; ++j)
      int_array_append((*tags)[tag_index], tag_index);
  }

  // Clean up.
  string_int_unordered_map_free(tag_indices);
  free(h_min);
}


