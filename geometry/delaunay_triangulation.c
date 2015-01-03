// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/array_utils.h"
#include "core/tuple.h"
#include "core/slist.h"
#include "core/unordered_set.h"
#include "geometry/tetrahedron.h"
#include "geometry/delaunay_triangulation.h"

// Algorithms for constructing Delaunay triangulations.
typedef enum
{
  BOWYER_WATSON,
  INCREMENTAL_FLIP,
  DIVIDE_AND_CONQUER,
  DEWALL
} delaunay_triangulation_algorithm_t;

// We use some of Shewchuk's robust geometric predicates.
extern real_t orient3d(real_t* pa, real_t* pb, real_t* pc, real_t* pd);
extern real_t insphere(real_t* pa, real_t* pb, real_t* pc, real_t* pd, real_t* pe);

struct delaunay_triangulation_t 
{
  delaunay_triangulation_algorithm_t algorithm;
  point_t* vertices;
  int num_vertices, vertex_cap, num_tets, tet_cap;
  int* tet_vertices;
};

// This helper allocates storage for a new vertex.
static void allocate_new_vertex(delaunay_triangulation_t* t)
{
  if ((t->num_vertices+1) >= t->vertex_cap)
  {
    while (t->vertex_cap < (t->num_vertices+1))
      t->vertex_cap *= 2;
    t->vertices = polymec_realloc(t->vertices, sizeof(point_t)*t->vertex_cap);
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
    t->tet_vertices = polymec_realloc(t->tet_vertices, 4*sizeof(int)*t->tet_cap);
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

// This helper performs a bistellar flip of two tetrahedra tau = (a, b, c, d) 
// and tau_a = (b, c, d, e) within the triangulation t. Notice that the 
// face common to tau and tau_a is (b, c, d).
static void flip(delaunay_triangulation_t* t, 
                 int tau, int tau_a, int a, int b, int c, int d, int e,
                 int_slist_t* stack)
{
  // Figure out how many facets of the tetrahedron tau_a are visible to 
  // point p within tetrahedron tau, and identify any degeneracies.
  point_t* p = &t->vertices[a];
  int tau_a_vertices[] = {b, c, d, e};
  static const int tau_a_facets[4][3] = {{1, 2, 3}, {0, 2, 3}, {0, 1, 3}, {0, 1, 2}};
  int num_visible_facets = 0;
  bool abdp_coplanar = false, p_lies_on_abc = false;
  for (int f = 0; f < 4; ++f)
  {
    point_t* v1 = &t->vertices[tau_a_vertices[tau_a_facets[f][0]]];
    point_t* v2 = &t->vertices[tau_a_vertices[tau_a_facets[f][1]]];
    point_t* v3 = &t->vertices[tau_a_vertices[tau_a_facets[f][2]]];
    // We use Shewchuk's ORIENT3D predicate to determine whether p lies above 
    // the face (v1, v2, v3).
    real_t pa[3] = {v1->x, v1->y, v1->z}, 
           pb[3] = {v2->x, v2->y, v2->z}, 
           pc[3] = {v3->x, v3->y, v3->z}, 
           pd[3] = {p->x, p->y, p->z};
    real_t orientation = orient3d(pa, pb, pc, pd);
    if (orientation < 0.0)
      ++num_visible_facets; // p is above (v1, v2, v3).
  }

  // Now flip tets according to the various cases in Ledoux's paper.
  if (num_visible_facets == 1) // case 1
  {
    // The union of tau and tau_a is convex, so a flip23 is performed.
    // FIXME
//    flip23(t, tau, tau_a);
//    int_slist_push(stack, ...);
  }
  else if (num_visible_facets == 2) // case 2
  {
    // The union of tau and tau_a is non-convex. Attempt to construct a 
    // tetrahedron tau_b = (a, b, p, d) such that the edge (a, b) is shared 
    // by tau, tau_a, and tau_b. If tau_b can be constructed, perform a flip32.
    bool abpd_exists = false;
    // FIXME
    if (abpd_exists)
    {
    }
    
    // If not, then no flip is performed, and this non-convex case will be 
    // rectified by another flip elsewhere.
  }
  else if (abdp_coplanar) // case 3
  {
    // a, b, d, and p are coplanar, so the creation of a tetrahedron (a, b, d, p) 
    // requires the tetrahedra tau and tau_a to be in the "config44" state.
    bool config44 = false;
    // FIXME
    if (config44)
    {
    }
  }
  else if (p_lies_on_abc) // case 4
  {
    // p lies on an edge of the face (a, b, c). So it is coplanar with both 
    // (a, b, c) and one other face. We perform a flip23 on tau and tau_a 
    // for reasons explained by Ledoux.
    // FIXME
  }
}

// This helper creates an initial set of tets within the triangulation so that 
// the convex hull covers the entire domain.
static void initial_aggregation(delaunay_triangulation_t* t, 
                                point_t* points, 
                                int num_points, 
                                int_unordered_set_t* points_in_t)
{
  // Find the points in the convex hull.
  // FIXME.
}


static void incremental_flip(delaunay_triangulation_t* t, point_t* points, int num_points)
{
  polymec_not_implemented("incremental_flip()");
  for (int i = 0; i < num_points; ++i)
  {
    // Insert the new vertex.
    allocate_new_vertex(t);
    int v_index = t->num_vertices;
    t->vertices[v_index] = points[i];
    t->num_vertices++;

    // Figure out which tet (tau) the point v lies in.
    int tau = tet_containing_point(t, &points[i]);

    // Insert a vertex into that tet and split it into 4 parts ("flip14" in 
    // Ledoux's paper).
    int new_tets[4];
    flip14(t, tau, new_tets);

    // Push the 4 new tetrahedra onto a working stack.
    int_slist_t* stack = int_slist_new();
    for (int j = 0; j < 4; ++j)
      int_slist_push(stack, new_tets[j]);

    // Perform all necessary flips to rectify the new vertex.
    while (!int_slist_empty(stack))
    {
      // Pull the tet tau off the stack and retrieve its vertex indices
      // (v, a, b, c), where v is the vertex we just inserted.
      int tau = int_slist_pop(stack, NULL);
      int v = -1, a = -1, b = -1, c = -1;
      for (int j = 0; j < 4; ++j)
      {
        if (t->tet_vertices[j] == v_index)
          v = t->tet_vertices[j];
        else
        {
          if (a == -1) 
            a = t->tet_vertices[j];
          else if (b == -1)
            b = t->tet_vertices[j];
          else if (c == -1)
            c = t->tet_vertices[j];
        }
      }

      // Find the tetrahedron (a, b, c, d) adjacent to tau = (v, a, b, c) 
      // that shares the vertices (a, b, c).
      int tau_a, d;
      find_adjacent_tet(t, tau, v, a, b, c, &tau_a, &d);

      // If the fourth vertex d of tau_a is inside the circumsphere of tau, 
      // flip tau and tau_a. Use Shewchuk's INSPHERE predicate for this.
      real_t pv[3] = {t->vertices[v].x, t->vertices[v].y, t->vertices[v].z};
      real_t pa[3] = {t->vertices[a].x, t->vertices[a].y, t->vertices[a].z};
      real_t pb[3] = {t->vertices[b].x, t->vertices[b].y, t->vertices[b].z};
      real_t pc[3] = {t->vertices[c].x, t->vertices[c].y, t->vertices[c].z};
      real_t pd[3] = {t->vertices[d].x, t->vertices[d].y, t->vertices[d].z};
      if (insphere(pv, pa, pb, pc, pd) > 0.0)
        flip(t, tau, tau_a, v, a, b, c, d, stack);
    }

    // Clean up.
    int_slist_free(stack);
  }
}

static void bowyer_watson(delaunay_triangulation_t* t, point_t* points, int num_points)
{
  // Initialize an aggregate of tets whose convex hull spans the entire domain.
  int_unordered_set_t* points_in_t = int_unordered_set_new();
  initial_aggregation(t, points, num_points, points_in_t);

  // Now add the rest of the points.
  tetrahedron_t* tet = tetrahedron_new();
  int_unordered_set_t* intersected_tets = int_unordered_set_new();
  int_tuple_unordered_set_t* intersected_faces = int_tuple_unordered_set_new();
  for (int i = 0; i < num_points; ++i)
  {
    // Skip vertices in the initial aggregation.
    if (int_unordered_set_contains(points_in_t, i))
      continue;

    point_t* x = &points[i];

    // Find all the tets whose circumcenters contain this point and 
    // mark them as intersected.
    for (int j = 0; j < t->num_tets; ++j)
    {
      point_t *xa = &points[t->tet_vertices[4*j]],
              *xb = &points[t->tet_vertices[4*j+1]],
              *xc = &points[t->tet_vertices[4*j+2]],
              *xd = &points[t->tet_vertices[4*j+3]];
      tetrahedron_set_vertices(tet, xa, xb, xc, xd);

      if (tetrahedron_contains_point(tet, x))
        int_unordered_set_insert(intersected_tets, j);
    }

    // Now make a list of the faces of these intersected tets. If a face 
    // occurs twice in this list, it is removed.
    int pos = 0, tet_index;
    while (int_unordered_set_next(intersected_tets, &pos, &tet_index))
    {
      // A face of a tet is stored as an ordered triple of its point indices.

      // This array allow us to express the 4 faces of this tet in terms 
      // of offsets from the first vertex of the tet.
      int offsets[4][3] = {{0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}};

      // Add each of these tet faces to our list. If it's already in there, 
      // take it out instead.
      for (int j = 0; j < 4; ++j)
      {
        int* f = int_tuple_new(3);
        f[0] = t->tet_vertices[4*tet_index+offsets[j][0]];
        f[1] = t->tet_vertices[4*tet_index+offsets[j][1]];
        f[2] = t->tet_vertices[4*tet_index+offsets[j][2]];
        int_qsort(f, 3);

        if (int_tuple_unordered_set_contains(intersected_faces, f))
          int_tuple_unordered_set_delete(intersected_faces, f);
        else
          int_tuple_unordered_set_insert(intersected_faces, f);
      }
    }

    // Now we assemble new tets consisting of the above faces.

    // First, overwrite any of the tets that were intersected.
    pos = 0;
    while (int_unordered_set_next(intersected_tets, &pos, &tet_index))
    {

    }

    // Now add new tets until we exhaust the rest of the intersected faces.
    while (!int_tuple_unordered_set_empty(intersected_faces))
    {
    }
  }

  // Clean up.
  int_unordered_set_free(intersected_tets);
  int_tuple_unordered_set_free(intersected_faces);
  int_unordered_set_free(points_in_t);
}

static void divide_and_conquer(delaunay_triangulation_t* t, point_t* points, int num_points)
{
  polymec_not_implemented("divide_and_conquer()");
}

static void dewall(delaunay_triangulation_t* t, point_t* points, int num_points)
{
  polymec_not_implemented("dewall()");
}

delaunay_triangulation_t* delaunay_triangulation_new(point_t* points, int num_points)
{
  ASSERT(num_points >= 4);

  delaunay_triangulation_t* t = polymec_malloc(sizeof(delaunay_triangulation_t));
  t->algorithm = BOWYER_WATSON;
  t->num_vertices = 0;
  t->vertex_cap = 32;
  t->vertices = polymec_malloc(sizeof(point_t) * t->vertex_cap);

  t->num_tets = 0;
  t->tet_cap = 32;
  t->tet_vertices = polymec_malloc(sizeof(int) * 4 * t->tet_cap);

  switch(t->algorithm)
  {
    case BOWYER_WATSON:
      bowyer_watson(t, points, num_points);
      break;
    case INCREMENTAL_FLIP:
      incremental_flip(t, points, num_points);
      break;
    case DIVIDE_AND_CONQUER:
      divide_and_conquer(t, points, num_points);
      break;
    case DEWALL:
      dewall(t, points, num_points);
  }

  return t;
}

void delaunay_triangulation_free(delaunay_triangulation_t* t)
{
  polymec_free(t->vertices);
  polymec_free(t->tet_vertices);
  polymec_free(t);
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

