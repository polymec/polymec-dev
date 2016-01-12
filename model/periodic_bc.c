// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <gc/gc.h>
#include "core/point.h"
#include "core/unordered_set.h"
#include "model/periodic_bc.h"

static const int periodic_bc_magic_number = 123652234;

// This function can be used by default to generate a periodic map.
static int_int_unordered_map_t* generate_periodic_map(void* context, mesh_t* mesh, char* tag1, char* tag2)
{
  // First, we validate the periodic mapping a bit.

  // Are there the same number of faces in both the periodic tags?
  int num_faces1, num_faces2;
  int* faces1 = mesh_tag(mesh->face_tags, (const char*)tag1, &num_faces1);
  int* faces2 = mesh_tag(mesh->face_tags, (const char*)tag2, &num_faces2);
  if (num_faces1 != num_faces2)
  {
    polymec_error("Number of faces for periodic boundary '%s' (%d)\n"
                  "does not equal number of faces for periodic boundary '%s' (%d)", num_faces1, tag1, num_faces2, tag2);
  }

  // All of the centers of the faces must be coplanar, so we pick three points 
  // in each and assemble a plane representation.
  point_t xp1, xp2;
  vector_t n1, n2;
  point_copy(&xp1, &mesh->face_centers[faces1[0]]);
  point_copy(&xp2, &mesh->face_centers[faces2[0]]);
  for (int p = 0; p < 2; ++p)
  {
    // Find 3 non-colinear points.
    point_t p1, p2, p3;
    point_copy(&p1, &mesh->face_centers[faces1[0]]);
    point_copy(&p2, &mesh->face_centers[faces1[1]]);
    point_copy(&p3, &mesh->face_centers[faces1[2]]);
    for (int i = 3; points_are_colinear(&p1, &p2, &p3); ++i)
      point_copy(&p3, &mesh->face_centers[faces1[i]]);

    // Find the normal vector of the plane containing these points.
    vector_t v1, v2, n;
    point_displacement(&p1, &p2, &v1);
    point_displacement(&p1, &p3, &v2);
    vector_cross(&v1, &v2, &n);
    vector_normalize(&n);
    if (p == 0)
      vector_copy(&n1, &n);
    else
      vector_copy(&n2, &n);
  }

  // Find the normal displacement vector D12 that maps a point from plane 1 
  // to plane 2.
  {
    vector_t D12;
    point_displacement(&xp1, &xp2, &D12);
    real_t D12_mag = vector_mag(&D12);
    real_t D12_n = vector_dot(&D12, &n1);
    if (D12_n < 0.0)
      vector_scale(&D12, -D12_n/D12_mag);
    else
      vector_scale(&D12, D12_n/D12_mag);
  }

  // Now that we are somewhat reassured of the sanity of the alleged 
  // periodicity, we can build the mapping. 
  int_unordered_set_t* face_mapped1 = int_unordered_set_new();
  int_unordered_set_t* face_mapped2 = int_unordered_set_new();
  int_int_unordered_map_t* pmap = int_int_unordered_map_new();
  for (int f = 0; f < num_faces1; ++f)
  {

  }
  // FIXME

  // Clean up.
  int_unordered_set_free(face_mapped1);
  int_unordered_set_free(face_mapped2);

  return pmap;
}

struct periodic_bc_t 
{
  // This magic number is used to validate periodic_bc_t objects which 
  // are cast from void pointers.
  int magic_number;

  // Essentially, this type contains the tags which are identified through 
  // a periodic boundary condition.
  char* tag1;
  char* tag2;

  // Function pointer (and context pointer) for generating periodic maps.
  periodic_bc_mapping_func generate_map;
  void* generate_map_context;
};

static void periodic_bc_free(void* ctx, void* dummy)
{
  periodic_bc_t* bc = (periodic_bc_t*)ctx;
  free(bc->tag1);
  free(bc->tag2);
}

periodic_bc_t* periodic_bc_new(const char* tag1, const char* tag2)
{
  // Use the default periodic mapping function.
  return periodic_bc_new_with_map_func(tag1, tag2, generate_periodic_map, NULL);
}

periodic_bc_t* periodic_bc_new_with_map_func(const char* tag1, const char* tag2, periodic_bc_mapping_func mapping_func, void* mapping_context)
{
  ASSERT(tag1 != NULL);
  ASSERT(tag2 != NULL);
  ASSERT(strcmp(tag1, tag2) != 0);

  periodic_bc_t* bc = GC_MALLOC(sizeof(periodic_bc_t));
  bc->magic_number = periodic_bc_magic_number;
  bc->tag1 = string_dup(tag1);
  bc->tag2 = string_dup(tag2);
  GC_register_finalizer(bc, &periodic_bc_free, bc, NULL, NULL);

  // Set up the map generation stuff.
  bc->generate_map = mapping_func;
  bc->generate_map_context = mapping_context;

  return bc;
}

bool periodic_bc_is_valid(periodic_bc_t* bc)
{
  return (bc->magic_number == periodic_bc_magic_number);
}

bool pointer_is_periodic_bc(void* ptr)
{
  periodic_bc_t* bc = (periodic_bc_t*)ptr;
  return periodic_bc_is_valid(bc);
}

void periodic_bc_get_tags(periodic_bc_t* bc, char** tag1, char** tag2)
{
  ASSERT(tag1 != NULL);
  ASSERT(tag2 != NULL);
  *tag1 = bc->tag1;
  *tag2 = bc->tag2;
}

int_int_unordered_map_t* periodic_bc_generate_map(periodic_bc_t* bc, mesh_t* mesh)
{
  // Just call the function pointer for this one.
  return bc->generate_map(bc->generate_map_context, mesh, bc->tag1, bc->tag2);
}



