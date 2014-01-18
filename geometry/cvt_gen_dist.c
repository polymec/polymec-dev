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

#include <gc/gc.h>
#include "geometry/cvt_gen_dist.h"
#include "geometry/scaled.h"

struct cvt_gen_dist_t 
{
  char* name;
  void* context;
  real_t scale_factor;
  cvt_gen_dist_vtable vtable;
};

static void cvt_gen_dist_free(void* ctx, void* dummy)
{
  cvt_gen_dist_t* dist = (cvt_gen_dist_t*)ctx;
  free(dist->name);
  if ((dist->context != NULL) && (dist->vtable.dtor != NULL))
    dist->vtable.dtor(dist->context);
}

static void project_point_to_boundary(sp_func_t* boundary, point_t* x)
{
  ASSERT(sp_func_has_deriv(boundary, 1));
  real_t D, grad_D[3];
  static int max_proj = 100;
  int i = 0;
  do 
  {
    sp_func_eval(boundary, x, &D);
    sp_func_eval_deriv(boundary, 1, x, grad_D);
    vector_t normal = {.x = -grad_D[0], .y = -grad_D[1], .z = -grad_D[2]};
    vector_normalize(&normal);
//printf("%d: %g %g %g (%g, %g, %g, %g) ->", i, x->x, x->y, x->z, D, normal.x, normal.y, normal.z);
    x->x += D * normal.x;
    x->y += D * normal.y;
    x->z += D * normal.z;
//printf("%g %g %g\n", x->x, x->y, x->z);
    ++i;
  }
  while ((D > 1e-12) && (i < max_proj));
  if (i == max_proj)
    polymec_error("project_to_boundary: Given boundary is not a signed distance function.");
}

cvt_gen_dist_t* cvt_gen_dist_new(const char* name, void* context, cvt_gen_dist_vtable vtable)
{
  ASSERT(vtable.iterate != NULL);

  cvt_gen_dist_t* dist = GC_MALLOC(sizeof(cvt_gen_dist_t));
  dist->name = string_dup(name);
  dist->context = context;
  dist->vtable = vtable;
  dist->scale_factor = 1.0;
  GC_register_finalizer(dist, cvt_gen_dist_free, dist, NULL, NULL);
  return dist;
}

const char* cvt_gen_dist_name(cvt_gen_dist_t* dist)
{
  return (const char*)dist->name;
}

static inline void swap_point(point_t* points, int i, int j)
{
  point_t temp;
  point_copy(&temp, &points[i]);
  point_copy(&points[i], &points[j]);
  point_copy(&points[j], &temp);
}

void cvt_gen_dist_iterate(cvt_gen_dist_t* dist, 
                          sp_func_t* density,
                          sp_func_t* boundary,
                          bbox_t* bounding_box,
                          point_t* points, 
                          int num_points,
                          int* num_boundary_points)
{
  ASSERT(num_points > 0);
  ASSERT(points != NULL);
  ASSERT(num_boundary_points != NULL);

  // If necessary, project points back to the boundary.
  if (boundary != NULL)
  {
    for (int i = 0; i < num_points; ++i)
    {
      real_t D;
      sp_func_eval(boundary, &points[i], &D);
      if (D >= 0.0)
      {
        project_point_to_boundary(boundary, &points[i]);
        real_t D;
        sp_func_eval(boundary, &points[i], &D);
        ASSERT(fabs(D) < 1e-12);
      }
    }
  }

  // Perform the iteration.
  dist->vtable.iterate(dist->context, density, boundary, bounding_box, 
                       points, num_points);

  // Check that all points are within the domain or on the boundary, and 
  // move the boundary points to the end of the list.
  *num_boundary_points = 0;
  if (boundary != NULL)
  {
    for (int i = 0; i < num_points - *num_boundary_points; ++i)
    {
      real_t D;
      sp_func_eval(boundary, &points[i], &D);
      if (D > 1e-12)
      {
        polymec_error("cvt_gen_dist_iterate: Point %d lies outside the domain (D = %g).", i, D);
      }
      else if (fabs(D) < 1e-12)
      {
        // This point is on the boundary, so we move it to the end of the list.
        ++(*num_boundary_points);
        swap_point(points, i, num_points - *num_boundary_points);
        --i;
      }
    }
  }
}

void cvt_gen_dist_generate_random_points(int (*rng)(), sp_func_t* density, bbox_t* bounding_box, int num_points, point_t* points)
{
  ASSERT(density != NULL);
  ASSERT(bounding_box != NULL);
  ASSERT(num_points > 0);
  ASSERT(points != NULL);

  // Keep track of the maximum value of the density we've hit so far.
  real_t rho_max = 1e-12;

  for (int i = 0; i < num_points; ++i)
  {
    bool rejected;
    do
    {
      // Generate a random point within the bounding box.
      point_randomize(&points[i], rng, bounding_box);

      // Now evaluate the density at this point and reject the point if 
      // if the relative density is less than a random number between 0 
      // and 1.
      real_t rho;
      sp_func_eval(density, &points[i], &rho);
      if (rho > rho_max) rho_max = rho;

      real_t P = rng()/RAND_MAX;
      rejected = (rho/rho_max < P);
    }
    while (rejected);
  }
}

