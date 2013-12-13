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


#include "core/unordered_set.h"
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
  iter->name = string_dup(name);
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
                   char** tag_names, int_array_t** tags, int num_tags,
                   cvt_iterator_t* cvt_iter)
{
  ASSERT(num_stationary_generators >= 0);
  ASSERT(num_mobile_generators > 0);
  ASSERT(((tag_names == NULL) && (tags == NULL) && (num_tags == 0)) || ((tag_names != NULL) && (tags != NULL) && (num_tags > 0)));
  ASSERT(cvt_iter != NULL);

  // Initialize our iterator if needed.
  if (cvt_iter->vtable.init != NULL)
  {
    cvt_iter->vtable.init(cvt_iter->context, stationary_generators, num_stationary_generators,
                          mobile_generators, num_mobile_generators);
  }

  // Create an initial tessellation from all the points.
  int num_generators = num_stationary_generators + num_mobile_generators;
  log_info("%s: creating initial tessellation (%d generators)", cvt_iter->name, num_generators);
  point_t* all_generators = malloc(sizeof(point_t) * num_generators);
  point_t* my_mobile_generators = malloc(sizeof(point_t) * num_mobile_generators);
  memcpy(all_generators, stationary_generators, sizeof(point_t) * num_stationary_generators);
  memcpy(&all_generators[num_stationary_generators], mobile_generators, sizeof(point_t) * num_mobile_generators);
  memcpy(my_mobile_generators, mobile_generators, sizeof(point_t) * num_mobile_generators);
  int_slist_t* deleted_generators = int_slist_new();
  mesh_t* mesh = create_voronoi_mesh(all_generators, num_generators, NULL, 0, deleted_generators);

  // Make sure that no mobile generators were deleted by the tessellation 
  // process.
  for (int_slist_node_t* dgen_iter = deleted_generators->front;
       dgen_iter != NULL; dgen_iter = dgen_iter->next)
  {
    if (dgen_iter->value > num_stationary_generators)
    {
      point_t* xg = &all_generators[dgen_iter->value];
      polymec_error("create_cvt: mobile generator %d at (%g, %g, %g) was deleted to bound the tessellation.\n"
                    "Please ensure that the mobile generators are bounded by stationary generators.", dgen_iter->value, xg->x, xg->y, xg->z);
    }
  }

  // Iterate till we're done.
  int iteration = 0;
  while (!cvt_iter->vtable.is_finished(cvt_iter->context, mesh, iteration))
  {
    log_info("%s: iteration %d", cvt_iter->name, iteration);

    // Move the mobile points.
    cvt_iter->vtable.move_points(cvt_iter->context, my_mobile_generators, num_mobile_generators);

    // Destroy and recreate the mesh.
    mesh_free(mesh);
    memcpy(all_generators, stationary_generators, sizeof(point_t) * num_stationary_generators);
    memcpy(&all_generators[num_stationary_generators], mobile_generators, sizeof(point_t) * num_mobile_generators);
    mesh = create_voronoi_mesh(all_generators, num_generators, NULL, 0, deleted_generators);

    // Make sure that no mobile generators were deleted by the tessellation 
    // process.
    for (int_slist_node_t* dgen_iter = deleted_generators->front;
         dgen_iter != NULL; dgen_iter = dgen_iter->next)
    {
      if (dgen_iter->value > num_stationary_generators)
      {
        polymec_error("create_cvt: mobile generator %d was deleted to bound the tessellation.\n"
                      "Please ensure that the mobile generators are bounded by stationary generators.");
      }
    }

    ++iteration;
  }

  // Clean up the stuff we're finished with.
  free(my_mobile_generators);
  free(all_generators);
  cvt_iterator_free(cvt_iter);

  // If we are given tags, add them to the mesh now. Note that we must 
  // account for any generators that were deleted during the tessellation.
  if (num_tags > 0)
  {
    // Create a mapping of generators to cell indices.
    int* cell_indices = malloc(sizeof(int) * num_generators);
    int_unordered_set_t* deleted_cells = int_unordered_set_new();
    for (int_slist_node_t* dgen_iter = deleted_generators->front;
         dgen_iter != NULL; dgen_iter = dgen_iter->next)
    {
      int_unordered_set_insert(deleted_cells, dgen_iter->value);
    }
    int cell_offset = 0;
    for (int i = 0; i < num_generators; ++i)
    {
      if (int_unordered_set_contains(deleted_cells, i))
        cell_indices[i] = -1; // Deleted generators are mapped to -1.
      else
      {
        cell_indices[i] = cell_offset;
        ++cell_offset;
      }
    }
    int_unordered_set_free(deleted_cells);

    // Now go through the tags and make a list of all the deleted cells in 
    // each one so that we can recompute the sizes.
    int* num_deleted_cells_in_tag = malloc(sizeof(int) * num_tags);
    memset(num_deleted_cells_in_tag, 0, sizeof(int) * num_tags);
    for (int i = 0; i < num_tags; ++i)
    {
      for (int j = 0; j < tags[i]->size; ++j)
      {
        ASSERT(tags[i]->data[j] >= 0);
        ASSERT(tags[i]->data[j] < num_generators);
        if (tags[i]->data[j] == -1)
          num_deleted_cells_in_tag[i] += 1;
      }
    }

    for (int i = 0; i < num_tags; ++i)
    {
      int tag_size = tags[i]->size - num_deleted_cells_in_tag[i];
      int* tagi = mesh_create_tag(mesh->cell_tags, tag_names[i], tag_size);
      int k = 0;
      for (int j = 0; j < tag_size; ++j)
      {
        while (cell_indices[tags[i]->data[j]] == -1) ++k;
        tagi[j] = cell_indices[k];
      }
    }
    free(num_deleted_cells_in_tag);
    free(cell_indices);
  }

  // Clean up the rest.
  int_slist_free(deleted_generators);

  return mesh;
}


