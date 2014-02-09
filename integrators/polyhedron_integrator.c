// Copyright (c) 2012-2014, Jeffrey N. Johnson
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

#include "geometry/polygon.h"
#include "integrators/polyhedron_integrator.h"

struct polyhedron_integrator_t 
{
  char* name;
  void* context;
  polyhedron_integrator_vtable vtable;
};

polyhedron_integrator_t* polyhedron_integrator_new(const char* name,
                                                   void* context,
                                                   polyhedron_integrator_vtable vtable)
{
  ASSERT(vtable.set_domain != NULL);
  ASSERT(vtable.next_volume_point != NULL);
  ASSERT(vtable.next_surface_point != NULL);
  polyhedron_integrator_t* integ = malloc(sizeof(polyhedron_integrator_t));
  integ->name = string_dup(name);
  integ->context = context;
  integ->vtable = vtable;
  return integ;
}

polyhedron_integrator_t* faceted_polyhedron_integrator_new(int order)
{
  return NULL;
}

polyhedron_integrator_t* cyl_wedge_1d_polyhedron_integrator_new(int order)
{
  return NULL;
}

polyhedron_integrator_t* cyl_1d_polyhedron_integrator_new(int order)
{
  return NULL;
}

polyhedron_integrator_t* sph_spike_1d_polyhedron_integrator_new(int order)
{
  return NULL;
}

polyhedron_integrator_t* sph_1d_polyhedron_integrator_new(int order)
{
  return NULL;
}

void polyhedron_integrator_free(polyhedron_integrator_t* integ)
{
  free(integ->name);
  if ((integ->vtable.dtor != NULL) && (integ->context != NULL))
    integ->vtable.dtor(integ->context);
  free(integ);
}

void polyhedron_integrator_set_domain(polyhedron_integrator_t* integ,
                                      point_t* face_nodes,
                                      int* face_node_offsets,
                                      int num_faces)
{
  integ->vtable.set_domain(integ->context, face_nodes, face_node_offsets, num_faces);
}

bool polyhedron_integrator_next_volume_point(polyhedron_integrator_t* integ,
                                             int* pos,
                                             point_t* point,
                                             real_t* weight)
{
  return integ->vtable.next_volume_point(integ->context, pos, point, weight);
}

bool polyhedron_integrator_next_surface_point(polyhedron_integrator_t* integ,
                                              int* pos,
                                              point_t* point,
                                              vector_t* normal_vector,
                                              real_t* weight)
{
  return integ->vtable.next_surface_point(integ->context, pos, point, normal_vector, weight);
}

typedef struct
{
  real_t cell_volume;
  point_t cell_center;

  int num_faces, face_cap;
  point_t* face_centers;
  real_t* face_areas;
  vector_t* face_normals;
} midpt_t;

static void midpt_set_domain(void* context, 
                             point_t* face_nodes, 
                             int* face_node_offsets,
                             int num_faces)
{
  midpt_t* midpt = context;
  midpt->num_faces = num_faces;
  if (midpt->num_faces > midpt->face_cap)
  {
    while (midpt->num_faces > midpt->face_cap)
      midpt->face_cap *= 2;
    midpt->face_centers = realloc(midpt->face_centers, sizeof(point_t) * midpt->face_cap);
    midpt->face_areas = realloc(midpt->face_areas, sizeof(real_t) * midpt->face_cap);
  }

  // Compute face centers/areas/normals and the cell center.
  polygon_t* polygons[midpt->num_faces];
  point_set(&midpt->cell_center, 0.0, 0.0, 0.0);
  for (int f = 0; f < midpt->num_faces; ++f)
  {
    int first_face_node = face_node_offsets[f];
    int num_face_nodes = face_node_offsets[f+1] - face_node_offsets[f];
    polygons[f] = polygon_new(&face_nodes[first_face_node], num_face_nodes);
    polygon_compute_centroid(polygons[f], &midpt->face_centers[f]);
    polygon_compute_normal(polygons[f], &midpt->face_normals[f]);
    midpt->face_areas[f] = polygon_area(polygons[f]);

    midpt->cell_center.x += midpt->face_centers[f].x;
    midpt->cell_center.y += midpt->face_centers[f].y;
    midpt->cell_center.z += midpt->face_centers[f].z;
  }
  midpt->cell_center.x /= midpt->num_faces;
  midpt->cell_center.y /= midpt->num_faces;
  midpt->cell_center.z /= midpt->num_faces;

  // Compute the polyhedron centroid and volume.
  // We compute the cell volume using the Divergence Theorem:
  // V = 1/3 * sum(f, xf o nf)
  midpt->cell_volume = 0.0;
  for (int f = 0; f < midpt->num_faces; ++f)
  {
    vector_t xf;
    point_displacement(&midpt->cell_center, &midpt->face_centers[f], &xf);
    midpt->cell_volume += 1.0/3.0 * vector_dot(&xf, &midpt->face_normals[f]);
  }
}

static bool midpt_next_volume_point(void* context,
                                    int* pos,
                                    point_t* point,
                                    real_t* weight)
{
  midpt_t* midpt = context;
  if (*pos == 0)
  {
    point_copy(point, &midpt->cell_center);
    *weight = midpt->cell_volume;
    ++(*pos);
    return true;
  }
  else
    return false;
}

static bool midpt_next_surface_point(void* context,
                                     int* pos,
                                     point_t* point,
                                     vector_t* normal_vector,
                                     real_t* weight)
{
  midpt_t* midpt = context;
  if (*pos < midpt->num_faces)
  {
    int f = *pos;
    point_copy(point, &midpt->face_centers[f]);
    vector_copy(normal_vector, &midpt->face_normals[f]);
    *weight = midpt->face_areas[f];
    ++(*pos);
    return true;
  }
  else
    return false;
}

static void midpt_dtor(void* context)
{
  midpt_t* midpt = context;
  free(midpt->face_normals);
  free(midpt->face_centers);
  free(midpt->face_areas);
  free(midpt);
}

polyhedron_integrator_t* midpoint_polyhedron_integrator_new()
{
  midpt_t* midpt = malloc(sizeof(midpt_t));
  midpt->cell_volume = 0.0;
  point_set(&midpt->cell_center, 0.0, 0.0, 0.0);
  midpt->face_cap = 12;
  midpt->face_centers = malloc(sizeof(point_t) * midpt->face_cap);
  memset(midpt->face_centers, 0, sizeof(point_t) * midpt->face_cap);
  midpt->face_areas = malloc(sizeof(real_t) * midpt->face_cap);
  memset(midpt->face_areas, 0, sizeof(real_t) * midpt->face_cap);
  midpt->face_normals = malloc(sizeof(vector_t) * midpt->face_cap);
  memset(midpt->face_normals, 0, sizeof(vector_t) * midpt->face_cap);

  polyhedron_integrator_vtable vtable = {.set_domain = midpt_set_domain,
                                         .next_volume_point = midpt_next_volume_point,
                                         .next_surface_point = midpt_next_surface_point,
                                         .dtor = midpt_dtor};
  return polyhedron_integrator_new("Midpoint", midpt, vtable);
}

