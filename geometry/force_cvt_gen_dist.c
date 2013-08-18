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

#include <gc/gc.h>
#include "core/slist.h"
#include "geometry/force_cvt_gen_dist.h"
#include "geometry/delaunay_triangulator.h"

// Types of forces.
typedef enum
{
  FORCE_LINEAR_SPRING
} cvt_gen_force_type_t;

struct cvt_gen_force_t 
{
  cvt_gen_force_type_t type;
  double param1;
};

cvt_gen_force_t* linear_spring_force_new(double k)
{
  cvt_gen_force_t* F = GC_MALLOC(sizeof(cvt_gen_force_t));
  F->type = FORCE_LINEAR_SPRING;
  F->param1 = k;
  return F;
}

typedef struct
{
  cvt_gen_force_t* force;
  double tolerance;
} force_cvt_gen_dist_t;

static void move_points(point_t* points, int num_points, vector_t* forces, double dt)
{
  for (int i = 0; i < num_points; ++i)
  {
    points[i].x += dt * forces[i].x;
    points[i].y += dt * forces[i].y;
    points[i].z += dt * forces[i].z;
  }
}

static void accumulate_forces_on_vertices(cvt_gen_force_t* force,
                                          sp_func_t* density, 
                                          delaunay_triangulation_t* t,
                                          double* forces,
                                          double* max_force)
{
  // FIXME: Assume linear spring model for now.
  ASSERT(force->type == FORCE_LINEAR_SPRING);
  double k = force->param1;

  // Compute equilibrium lengths for each of the vertices.
  point_t vertices[t->num_vertices];
  double Leq[t->num_vertices];
  for (int i = 0; i < t->num_vertices; ++i)
  {
    vertices[i].x = t->vertices[3*i];
    vertices[i].y = t->vertices[3*i+1];
    vertices[i].z = t->vertices[3*i+2];
    sp_func_eval(density, &vertices[i], &Leq[i]);
    Leq[i] = 1.0/sqrt(Leq[i]);
  }

  *max_force = 0.0;
  for (int i = 0; i < t->num_tets; ++i)
  {
    int j0 = t->tets[i].vertices[0];
    int j1 = t->tets[i].vertices[1];
    int j2 = t->tets[i].vertices[2];
    int j3 = t->tets[i].vertices[3];

    // Distances.
    double L01 = point_distance(&vertices[j0], &vertices[j1]);
    double L02 = point_distance(&vertices[j0], &vertices[j2]);
    double L03 = point_distance(&vertices[j0], &vertices[j3]);
    double L12 = point_distance(&vertices[j1], &vertices[j2]);
    double L13 = point_distance(&vertices[j1], &vertices[j3]);
    double L23 = point_distance(&vertices[j2], &vertices[j3]);

    // Force magnitudes.
    double L0 = Leq[j0];
    double F01 = MAX(k * (L0 - L01), 0.0);
    double F02 = MAX(k * (L0 - L02), 0.0);
    double F03 = MAX(k * (L0 - L03), 0.0);

    double L1 = Leq[j1];
    double F10 = MAX(k * (L1 - L01), 0.0);
    double F12 = MAX(k * (L1 - L12), 0.0);
    double F13 = MAX(k * (L1 - L13), 0.0);

    double L2 = Leq[j2];
    double F20 = MAX(k * (L2 - L02), 0.0);
    double F21 = MAX(k * (L2 - L12), 0.0);
    double F23 = MAX(k * (L2 - L23), 0.0);

    double L3 = Leq[j3];
    double F30 = MAX(k * (L3 - L03), 0.0);
    double F31 = MAX(k * (L3 - L13), 0.0);
    double F33 = MAX(k * (L3 - L23), 0.0);
  }
}

void force_cvt_gen_dist_iterate(void* context, 
                                sp_func_t* density,
                                sp_func_t* boundary,
                                bbox_t* bounding_box,
                                point_t* points, 
                                int num_points)
{
  ASSERT(density != NULL);
  ASSERT(bounding_box != NULL);
  ASSERT(num_points > 0);

  force_cvt_gen_dist_t* algo = context;

  delaunay_triangulator_t* triangulator = delaunay_triangulator_new();

  vector_t forces[num_points];

  // Create an initial triangulation for these points.
  delaunay_triangulation_t* t = delaunay_triangulator_triangulate(triangulator, points, num_points);

  // Compute the initial forces on each point.
  for (int i = 0; i < num_points; ++i)
    forces[i].x = forces[i].y = forces[i].z = 0.0;
  double Fmax;
  accumulate_forces_on_vertices(algo->force, density, t, forces, &Fmax);
  double dt = 0.1;

  while (Fmax > algo->tolerance)
  {
    // Move the points according to the forces.
    move_points(points, num_points, forces, dt);

    // Recreate the triangulation for these points.
    delaunay_triangulation_free(t);
    t = delaunay_triangulator_triangulate(triangulator, points, num_points);

    // Now recompute the forces on each point.
    for (int i = 0; i < num_points; ++i)
      forces[i].x = forces[i].y = forces[i].z = 0.0;
    accumulate_forces_on_vertices(algo->force, density, t, forces, &Fmax);
  }
  while (Fmax > algo->tolerance);
}

cvt_gen_dist_t* force_cvt_gen_dist_new(cvt_gen_force_t* force,
                                       double tolerance)
{
  ASSERT(force != NULL);
  ASSERT(tolerance > 0.0);

  force_cvt_gen_dist_t* algo = malloc(sizeof(force_cvt_gen_dist_t));
  algo->force = force;
  algo->tolerance = tolerance;
  cvt_gen_dist_vtable vtable = {.iterate = force_cvt_gen_dist_iterate, .dtor = free};
  return cvt_gen_dist_new("Force balance CVT generator distribution", algo, vtable);
}

