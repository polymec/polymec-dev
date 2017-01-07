// Copyright (c) 2012-2017, Jeffrey N. Johnson
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
  int num_vertices;

  polygon_t** faces;
  int num_faces;
};

static void polyhedron_free(void* ctx)
{
  polyhedron_t* poly = ctx;
  polymec_free(poly->vertices);
  for (int i = 0; i < poly->num_faces; ++i)
    poly->faces[i] = NULL;
}

polyhedron_t* polyhedron_new(point_t* vertices, int num_vertices,
                             int** faces, int* num_face_vertices, 
                             int num_faces)
{
  polyhedron_t* poly = polymec_gc_malloc(sizeof(polyhedron_t), polyhedron_free);

  // Copy in the vertices.
  poly->num_vertices = num_vertices;
  poly->vertices = polymec_malloc(sizeof(point_t) * num_vertices);
  memcpy(poly->vertices, vertices, sizeof(point_t) * num_vertices);

  // Now construct polygonal faces from the vertices.
  poly->num_faces = num_faces;
  poly->faces = polymec_malloc(sizeof(polygon_t*) * num_faces);
  for (int f = 0; f < num_faces; ++f)
  {
    point_t face_vertices[num_face_vertices[f]];
    for (int v = 0; v < num_face_vertices[f]; ++v)
      face_vertices[v] = vertices[faces[f][v]];
    poly->faces[f] = polygon_new(face_vertices, num_face_vertices[f]);
  }

  return poly;
}

int polyhedron_num_vertices(polyhedron_t* poly)
{
  return poly->num_vertices;
}

int polyhedron_num_faces(polyhedron_t* poly)
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

bool polyhedron_next_face(polyhedron_t* poly, int* pos, polygon_t** face)
{
  if (*pos < poly->num_faces)
  {
    *face = poly->faces[*pos];
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
  for (int f = 0; f < poly->num_faces; ++f)
  {
    int pos = 0;
    point_t* v;
    polygon_t* face = poly->faces[f];
    int nv = polygon_num_vertices(face);
    point_t vertices[nv];
    while (polygon_next_vertex(face, &pos, &v))
      vertices[pos-1] = *v;
    point_t xf;
    polygon_compute_centroid(face, &xf);

    // This face has nv tets whose volumes we sum.
    for (int i = 0; i < nv; ++i)
    {
      point_t p1 = vertices[i];
      point_t p2 = vertices[(i+1)%nv];
      tetrahedron_set_vertices(tet, &p1, &p2, &xf, &xc);
      V += tetrahedron_volume(tet);
    }
  }
  return V;
}

void polyhedron_compute_centroid(polyhedron_t* poly, point_t* centroid)
{
  centroid->x = centroid->y = centroid->z = 0.0;
  for (int f = 0; f < poly->num_faces; ++f)
  {
    point_t xf;
    polygon_compute_centroid(poly->faces[f], &xf);
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
  clone->vertices = polymec_malloc(sizeof(point_t*) * poly->num_vertices);
  memcpy(clone->vertices, poly->vertices, sizeof(point_t*) * poly->num_vertices);
  clone->num_faces = poly->num_faces;
  clone->faces = polymec_malloc(sizeof(polygon_t*) * poly->num_faces);
  memcpy(clone->faces, poly->faces, sizeof(polygon_t*) * poly->num_faces);
  return clone;
}

