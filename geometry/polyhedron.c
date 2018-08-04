// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "geometry/polyhedron.h"
#include "geometry/tetrahedron.h"

struct polyhedron_t 
{
  point_t* vertices;
  size_t num_vertices;
  point_t** face_vertices;
  size_t* num_face_vertices;
  size_t num_faces;
};

static void polyhedron_free(void* ctx)
{
  polyhedron_t* poly = ctx;
  for (size_t i = 0; i < poly->num_faces; ++i)
    polymec_free(poly->face_vertices[i]);
  polymec_free(poly->num_face_vertices);
  polymec_free(poly->face_vertices);
  polymec_free(poly->vertices);
}

polyhedron_t* polyhedron_new(point_t* vertices, size_t num_vertices,
                             int** faces, size_t* num_face_vertices, 
                             size_t num_faces)
{
  polyhedron_t* poly = polymec_gc_malloc(sizeof(polyhedron_t), polyhedron_free);

  // Copy in the vertices.
  poly->num_vertices = num_vertices;
  poly->vertices = polymec_malloc(sizeof(point_t) * num_vertices);
  memcpy(poly->vertices, vertices, sizeof(point_t) * num_vertices);

  // Copy in face vertex info.
  poly->num_faces = num_faces;
  poly->face_vertices = polymec_malloc(sizeof(point_t*) * num_faces);
  poly->num_face_vertices = polymec_malloc(sizeof(size_t) * num_faces);
  memcpy(poly->num_face_vertices, num_face_vertices, sizeof(size_t) * num_faces);
  for (size_t f = 0; f < num_faces; ++f)
  {
    poly->face_vertices[f] = polymec_malloc(sizeof(point_t) * num_face_vertices[f]);
    for (size_t v = 0; v < num_face_vertices[f]; ++v)
      poly->face_vertices[f][v] = vertices[faces[f][v]];
  }

  return poly;
}

size_t polyhedron_num_vertices(polyhedron_t* poly)
{
  return poly->num_vertices;
}

size_t polyhedron_num_faces(polyhedron_t* poly)
{
  return poly->num_faces;
}

bool polyhedron_next_vertex(polyhedron_t* poly, int* pos, point_t* vertex)
{
  if (*pos < poly->num_vertices)
  {
    *vertex = poly->vertices[*pos];
    ++(*pos);
    return true;
  }
  else
    return false;
}

bool polyhedron_next_face(polyhedron_t* poly, int* pos, 
                          point_t** face_vertices, size_t* num_face_vertices)
{
  if (*pos < (int)poly->num_faces)
  {
    *face_vertices = poly->face_vertices;
    *num_face_vertices = poly->num_face_vertices;
    ++(*pos);
    return true;
  }
  else
    return false;
}

real_t polyhedron_volume(polyhedron_t* poly)
{
  real_t V = 0.0;

  // Sum up the volumes of the tets whose bases are triangulations 
  // of our faces, with tips at our centroid.
  point_t xc;
  polyhedron_compute_centroid(poly, &xc);

  tetrahedron_t* tet = tetrahedron_new();
  for (size_t f = 0; f < poly->num_faces; ++f)
  {
    // Compute the (3D) centroid of this face. 
    // FIXME
    point_t xf;
    polygon_compute_centroid(face, &xf);

    // Sum the volumes of the tets formed by this face's vertex and its 
    // centroid.
    size_t nv = poly->num_face_vertices[f];
    for (size_t i = 0; i < nv; ++i)
    {
      point_t p1 = poly->face_vertices[f][i];
      point_t p2 = poly->face_vertices[f][(i+1)%nv];
      tetrahedron_set_vertices(tet, &p1, &p2, &xf, &xc);
      V += tetrahedron_volume(tet);
    }
  }
  return V;
}

void polyhedron_compute_centroid(polyhedron_t* poly, point_t* centroid)
{
  centroid->x = centroid->y = centroid->z = 0.0;
  for (size_t f = 0; f < poly->num_faces; ++f)
  {
    point_t xf;
    polygon_compute_centroid(poly->face_vertices[f], &xf);
    centroid->x += xf.x;
    centroid->y += xf.y;
    centroid->z += xf.z;
  }
  centroid->x /= poly->num_faces;
  centroid->y /= poly->num_faces;
  centroid->z /= poly->num_faces;
}

polyhedron_t* polyhedron_clone(polyhedron_t* poly)
{
  polyhedron_t* clone = polymec_gc_malloc(sizeof(polyhedron_t), polyhedron_free);
  clone->num_vertices = poly->num_vertices;
  clone->vertices = polymec_malloc(sizeof(point_t) * poly->num_vertices);
  memcpy(clone->vertices, poly->vertices, sizeof(point_t) * poly->num_vertices);

  clone->num_faces = poly->num_faces;
  clone->face_vertices = polymec_malloc(sizeof(point_t*) * poly->num_faces);
  clone->num_face_vertices = polymec_malloc(sizeof(size_t) * poly->num_faces);
  memcpy(clone->num_face_vertices, poly->num_face_vertices, sizeof(size_t) * poly->num_faces);
  for (size_t f = 0; f < poly->num_faces; ++f)
  {
    clone->face_vertices[f] = polymec_malloc(sizeof(point_t) * poly->num_face_vertices[f]);
    memcpy(clone->face_vertices[f], poly->face_vertices[f], sizeof(point_t) * poly->num_face_vertices[f]);
  }
  return clone;
}

