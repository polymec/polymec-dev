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

#include "core/slist.h"
#include "geometry/tetrahedron.h"
#include "geometry/delaunay_triangulation.h"

struct delaunay_triangulation_t 
{
  point_t* vertices;
  int num_vertices, vertex_cap, num_tets, tet_cap;
  int* tet_vertices;
  tetrahedron_t* big_tet;
};

delaunay_triangulation_t* delaunay_triangulation_new(point_t* v1, point_t* v2, point_t* v3, point_t* v4)
{
  delaunay_triangulation_t* t = malloc(sizeof(delaunay_triangulation_t));
  t->num_vertices = 4;
  t->vertex_cap = 32;
  t->vertices = malloc(sizeof(point_t) * t->vertex_cap);
  t->vertices[0] = *v1;
  t->vertices[1] = *v2;
  t->vertices[2] = *v3;
  t->vertices[3] = *v4;

  t->num_tets = 1;
  t->tet_cap = 32;
  t->tet_vertices = malloc(sizeof(int) * 4 * t->tet_cap);
  t->tet_vertices[0] = 0;
  t->tet_vertices[1] = 1;
  t->tet_vertices[2] = 2;
  t->tet_vertices[3] = 3;

  // Set up the Big Tet.
  t->big_tet = tetrahedron_new();
  tetrahedron_set_vertices(t->big_tet, v1, v2, v3, v4);

  return t;
}

void delaunay_triangulation_free(delaunay_triangulation_t* t)
{
  t->big_tet = NULL;
  free(t->vertices);
  free(t->tet_vertices);
  free(t);
}

// This helper allocates storage for a new vertex.
static void allocate_new_vertex(delaunay_triangulation_t* t)
{
  if ((t->num_vertices+1) >= t->vertex_cap)
  {
    while (t->vertex_cap < (t->num_vertices+1))
      t->vertex_cap *= 2;
    t->vertices = realloc(t->vertices, sizeof(point_t)*t->vertex_cap);
  }
}

// This helper allocates storage for the given number of new tets.
static void allocate_new_tets(delaunay_triangulation_t* t, int num_new_tets)
{
  ASSERT(num_new_tets > 0);

  if ((t->num_tets+num_new_tets) >= t->tet_cap)
  {
    while (t->tet_cap < (t->num_tets+num_new_tets))
      t->tet_cap *= 2;
    t->tet_vertices = realloc(t->tet_vertices, 4*sizeof(int)*t->tet_cap);
  }
}

// This helper returns the index of the tetrahedron containing the point p.
static int tet_containing_point(delaunay_triangulation_t* t, point_t* p)
{
  // FIXME
  return 0;
}

// This helper inserts a vertex into the tetrahedron with the given index
// in the triangulation, breaking it into 4 tetrahedra. The indices of the 
// 4 new tetrahedra are stored in new_tet_indices.
static void flip14(delaunay_triangulation_t* t, int old_tet_index, int* new_tet_indices)
{
}

// Given a tetrahedron tau1 = (v, a, b, c), this helper retrieves the index 
// and the vertex d of the tetrahedron tau2 = (a, b, c, d) adjacent to tau1 
// that shares the vertices (a, b, c).
static void find_adjacent_tet(delaunay_triangulation_t* t, 
                              int tau1, int v, int a, int b, int c,
                              int* tau2, int* d)
{
}

// This helper performs a bistellar flip of two tetrahedra tau and tau_a within 
// the triangulation t.
static void flip(delaunay_triangulation_t* t, 
                 int tau, int tau_a, int v, int a, int b, int c, int d,
                 int_slist_t* stack)
{
  // Figure out how many facets of the tetrahedron tau_a are visible to 
  // point p within tetrahedron tau.
  point_t* p = &t->vertices[v];
  int tau_a_vertices[] = {a, b, c, d};
  static const int tau_a_facets[4][3] = {{1, 2, 3}, {0, 2, 3}, {0, 1, 3}, {0, 1, 2}};
  int num_visible_facets = 0;
  for (int f = 0; f < 4; ++f)
  {
    point_t* v1 = &t->vertices[tau_a_vertices[tau_a_facets[f][0]]];
    point_t* v2 = &t->vertices[tau_a_vertices[tau_a_facets[f][1]]];
    point_t* v3 = &t->vertices[tau_a_vertices[tau_a_facets[f][2]]];
  }

  if (num_visible_facets == 1) // case 1
  {
//    flip23(t, tau, tau_a);
//    int_slist_push(stack, ...);
  }
}

void delaunay_triangulation_insert_vertex(delaunay_triangulation_t* t, 
                                          point_t* v)
{
  ASSERT(tetrahedron_contains_point(t->big_tet, v));

  // Insert the new vertex.
  allocate_new_vertex(t);
  int v_index = t->num_vertices;
  t->vertices[v_index] = *v;
  t->num_vertices++;

  // Figure out which tet (tau) the point v lies in.
  int tau = tet_containing_point(t, v);

  // Insert a vertex into that tet and split it into 4 parts ("flip14" in 
  // Ledoux's paper).
  int new_tets[4];
  flip14(t, tau, new_tets);

  // Push the 4 new tetrahedra onto a working stack.
  int_slist_t* stack = int_slist_new();
  for (int i = 0; i < 4; ++i)
    int_slist_push(stack, new_tets[i]);

  // Perform all necessary flips to rectify the new vertex.
  while (!int_slist_empty(stack))
  {
    // Pull the tet tau off the stack and retrieve its vertex indices
    // (v, a, b, c), where v is the vertex we just inserted.
    int tau = int_slist_pop(stack, NULL);
    int v = -1, a = -1, b = -1, c = -1;
    for (int i = 0; i < 4; ++i)
    {
      if (t->tet_vertices[i] == v_index)
        v = t->tet_vertices[i];
      else
      {
        if (a == -1) 
          a = t->tet_vertices[i];
        else if (b == -1)
          b = t->tet_vertices[i];
        else if (c == -1)
          c = t->tet_vertices[i];
      }
    }

    // Find the tetrahedron (a, b, c, d) adjacent to tau = (v, a, b, c) 
    // that shares the vertices (a, b, c).
    int tau_a, d;
    find_adjacent_tet(t, tau, v, a, b, c, &tau_a, &d);

    // Compute the radius of the circumsphere of the tet tau.
    point_t tau_nodes[4], tau_xc;
    int tau_vertices[] = {v, a, b, c};
    delaunay_triangulation_get_vertices(t, tau_vertices, 4, tau_nodes);
    tetrahedron_t* tau_tet = tetrahedron_new();
    tetrahedron_set_vertices(tau_tet, &tau_nodes[0], &tau_nodes[1], &tau_nodes[2], &tau_nodes[3]);
    tetrahedron_compute_circumcenter(tau_tet, &tau_xc);
    real_t tau_r = point_distance(&tau_xc, &tau_nodes[0]); // radius of circumsphere.

    // If the fourth vertex of tau_a is inside the circumsphere of tau, 
    // flip tau and tau_a.
    point_t xd = t->vertices[d];
    if (point_distance(&tau_xc, &xd) < tau_r)
      flip(t, tau, tau_a, v, a, b, c, d, stack);
  }
}

int delaunay_triangulation_num_vertices(delaunay_triangulation_t* t)
{
  return t->num_vertices;
}

void delaunay_triangulation_get_vertices(delaunay_triangulation_t* t, int* indices, int num_vertices, point_t* vertices)
{
#ifndef NDEBUG
  for (int i = 0; i < num_vertices; ++i)
  {
    ASSERT(indices[i] >= 0);
    ASSERT(indices[i] < t->num_vertices);
  }
#endif
  for (int i = 0; i < num_vertices; ++i)
    vertices[i] = t->vertices[indices[i]];
}

int delaunay_triangulation_num_tetrahedra(delaunay_triangulation_t* t)
{
  return t->num_tets;
}

bool delaunay_triangulation_next(delaunay_triangulation_t* t,
                                 int* pos, int* v1, int* v2, int* v3, int* v4)
{
  if (*pos >= t->num_tets)
    return false;
  int i = *pos;
  *v1 = t->tet_vertices[4*i];
  *v2 = t->tet_vertices[4*i+1];
  *v3 = t->tet_vertices[4*i+2];
  *v4 = t->tet_vertices[4*i+3];
  ++(*pos);
  return true;
}

