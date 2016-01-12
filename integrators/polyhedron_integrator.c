// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

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
  polyhedron_integrator_t* integ = polymec_malloc(sizeof(polyhedron_integrator_t));
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
  polymec_free(integ->name);
  if ((integ->vtable.dtor != NULL) && (integ->context != NULL))
    integ->vtable.dtor(integ->context);
  polymec_free(integ);
}

void polyhedron_integrator_set_domain(polyhedron_integrator_t* integ,
                                      mesh_t* mesh,
                                      int cell)
{
  integ->vtable.set_domain(integ->context, mesh, cell);
}

bool polyhedron_integrator_next_volume_point(polyhedron_integrator_t* integ,
                                             int* pos,
                                             point_t* point,
                                             real_t* weight)
{
  return integ->vtable.next_volume_point(integ->context, pos, point, weight);
}

int polyhedron_integrator_num_volume_points(polyhedron_integrator_t* integ)
{
  int pos = 0, num_points = 0;
  point_t x;
  real_t w;
  while (polyhedron_integrator_next_volume_point(integ, &pos, &x, &w))
    ++num_points;
  return num_points;
}

bool polyhedron_integrator_next_surface_point(polyhedron_integrator_t* integ,
                                              int face,
                                              int* pos,
                                              point_t* point,
                                              vector_t* normal_vector,
                                              real_t* weight)
{
  return integ->vtable.next_surface_point(integ->context, face, pos, point, normal_vector, weight);
}

int polyhedron_integrator_num_surface_points(polyhedron_integrator_t* integ,
                                             int face)
{
  int pos = 0, num_points = 0;
  point_t x;
  vector_t n;
  real_t w;
  while (polyhedron_integrator_next_surface_point(integ, face, &pos, &x, &n, &w))
    ++num_points;
  return num_points;
}

typedef struct
{
  real_t cell_volume;
  point_t cell_center;

  int num_faces;
  point_t face_centers[50];
  real_t face_areas[50];
  vector_t face_normals[50];
} midpt_t;

static void midpt_set_domain(void* context, 
                             mesh_t* mesh, 
                             int cell)
{
  midpt_t* midpt = context;
  ASSERT(mesh_cell_num_faces(mesh, cell) < 50);
  midpt->num_faces = mesh_cell_num_faces(mesh, cell);

  // Fetch face centers/areas/normals and the cell center.
  midpt->cell_center = mesh->cell_centers[cell];
  midpt->cell_volume = mesh->cell_volumes[cell];
  int pos = 0, face;
  point_t* xc = &midpt->cell_center;
  while (mesh_cell_next_oriented_face(mesh, cell, &pos, &face))
  {
    int local_face = pos - 1;
    int actual_face = (face >= 0) ? face : ~face;
    midpt->face_areas[local_face] = mesh->face_areas[actual_face];
    midpt->face_centers[local_face] = mesh->face_centers[actual_face];

    // Use the orientation of the face to determine the direction of the 
    // face's normal vector.
    vector_t n = mesh->face_normals[actual_face];
    point_t* xf = &mesh->face_centers[actual_face];
    vector_t outward;
    point_displacement(xc, xf, &outward);
    if (face != actual_face)
      vector_scale(&n, -1.0);
    midpt->face_normals[local_face] = n;
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
                                     int face,
                                     int* pos,
                                     point_t* point,
                                     vector_t* normal_vector,
                                     real_t* weight)
{
  midpt_t* midpt = context;
  if (*pos == 0)
  {
    point_copy(point, &midpt->face_centers[face]);
    vector_copy(normal_vector, &midpt->face_normals[face]);
    *weight = midpt->face_areas[face];
    ++(*pos);
    return true;
  }
  else
    return false;
}

static void midpt_dtor(void* context)
{
  midpt_t* midpt = context;
  polymec_free(midpt);
}

polyhedron_integrator_t* midpoint_polyhedron_integrator_new()
{
  midpt_t* midpt = polymec_malloc(sizeof(midpt_t));
  midpt->cell_volume = 0.0;
  point_set(&midpt->cell_center, 0.0, 0.0, 0.0);

  polyhedron_integrator_vtable vtable = {.set_domain = midpt_set_domain,
                                         .next_volume_point = midpt_next_volume_point,
                                         .next_surface_point = midpt_next_surface_point,
                                         .dtor = midpt_dtor};
  return polyhedron_integrator_new("Midpoint", midpt, vtable);
}

