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
  real_t param1;
};

cvt_gen_force_t* linear_spring_force_new(real_t k)
{
  cvt_gen_force_t* F = GC_MALLOC(sizeof(cvt_gen_force_t));
  F->type = FORCE_LINEAR_SPRING;
  F->param1 = k;
  return F;
}

typedef struct
{
  cvt_gen_force_t* force;
  real_t tolerance;
} force_cvt_gen_dist_t;

static void move_points(point_t* points, int num_points, vector_t* forces, real_t dt)
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
                                          real_t* forces,
                                          real_t* max_force)
{
  // FIXME: Assume linear spring model for now.
  ASSERT(force->type == FORCE_LINEAR_SPRING);
  real_t k = force->param1;

  // Compute equilibrium lengths for each of the vertices.
  point_t vertices[t->num_vertices];
  real_t Leq[t->num_vertices];
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
    real_t L01 = point_distance(&vertices[j0], &vertices[j1]);
    real_t L02 = point_distance(&vertices[j0], &vertices[j2]);
    real_t L03 = point_distance(&vertices[j0], &vertices[j3]);
    real_t L12 = point_distance(&vertices[j1], &vertices[j2]);
    real_t L13 = point_distance(&vertices[j1], &vertices[j3]);
    real_t L23 = point_distance(&vertices[j2], &vertices[j3]);

    // Force magnitudes.
    real_t L0 = Leq[j0];
    real_t F01 = MAX(k * (L0 - L01), 0.0);
    real_t F02 = MAX(k * (L0 - L02), 0.0);
    real_t F03 = MAX(k * (L0 - L03), 0.0);

    real_t L1 = Leq[j1];
    real_t F10 = MAX(k * (L1 - L01), 0.0);
    real_t F12 = MAX(k * (L1 - L12), 0.0);
    real_t F13 = MAX(k * (L1 - L13), 0.0);

    real_t L2 = Leq[j2];
    real_t F20 = MAX(k * (L2 - L02), 0.0);
    real_t F21 = MAX(k * (L2 - L12), 0.0);
    real_t F23 = MAX(k * (L2 - L23), 0.0);

    real_t L3 = Leq[j3];
    real_t F30 = MAX(k * (L3 - L03), 0.0);
    real_t F31 = MAX(k * (L3 - L13), 0.0);
    real_t F33 = MAX(k * (L3 - L23), 0.0);
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
  real_t Fmax;
  accumulate_forces_on_vertices(algo->force, density, t, forces, &Fmax);
  real_t dt = 0.1;

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
                                       real_t tolerance)
{
  ASSERT(force != NULL);
  ASSERT(tolerance > 0.0);

  force_cvt_gen_dist_t* algo = malloc(sizeof(force_cvt_gen_dist_t));
  algo->force = force;
  algo->tolerance = tolerance;
  cvt_gen_dist_vtable vtable = {.iterate = force_cvt_gen_dist_iterate, .dtor = free};
  return cvt_gen_dist_new("Force balance CVT generator distribution", algo, vtable);
}

