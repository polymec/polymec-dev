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
#include "core/kd_tree.h"
#include "geometry/prob_cvt_gen_dist.h"
#include "geometry/generate_random_points.h"

typedef struct 
{
  int num_samples; // Number of sample points.
  real_t alpha1, beta1, alpha2, beta2; // Algorithm coefficients.
  int (*rng)(); // Random number generator.
  real_t min_dist; // Minimum distance of a point from the boundary.
  int max_iters; // Maximum number of iterations.
} prob_cvt_gen_dist_t;

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

static void correct_generator(real_t alpha1, real_t beta1, 
                              real_t alpha2, real_t beta2, 
                              ptr_slist_t* nearest,
                              point_t* zi, int* ji_ptr)
{
  if (nearest->size > 0)
  {
    // Compute the average, ui, of the sample points in the Voronoi region.
    point_t ui = {.x = 0.0, .y = 0.0, .z = 0.0};
    for (ptr_slist_node_t* node = nearest->front; node != NULL; node = node->next)
    {
      point_t* pos = node->value;
      ui.x += pos->x;
      ui.y += pos->y;
      ui.z += pos->z;
    }
    ui.x /= nearest->size;
    ui.y /= nearest->size;
    ui.z /= nearest->size;

    // Correct the position of the generator.
    int ji = *ji_ptr;
    zi->x = ((alpha1*ji + beta1)*zi->x + (alpha2*ji + beta2)*ui.x) / (ji + 1.0);
    zi->y = ((alpha1*ji + beta1)*zi->y + (alpha2*ji + beta2)*ui.y) / (ji + 1.0);
    zi->z = ((alpha1*ji + beta1)*zi->z + (alpha2*ji + beta2)*ui.z) / (ji + 1.0);

    // Increment ji.
    ++(*ji_ptr);
  }
}

void prob_cvt_gen_dist_iterate(void* context, 
                               sp_func_t* density,
                               sp_func_t* boundary,
                               bbox_t* bounding_box,
                               point_t* points, 
                               int num_points)
{
  ASSERT(density != NULL);
  ASSERT(bounding_box != NULL);
  ASSERT(num_points > 0);

  prob_cvt_gen_dist_t* prob = context;

  int iter = 0;

  // Set ji to 1 for all i.
  int ji[num_points];
  for (int i = 0; i < num_points; ++i)
    ji[i] = 1;

  real_t alpha1 = prob->alpha1;
  real_t beta1 = prob->beta1;
  real_t alpha2 = prob->alpha2;
  real_t beta2 = prob->beta2;

  // Iterate until termination.
  point_t samples[prob->num_samples];
  ptr_slist_t* near_ipoints[num_points];
  for (int i = 0; i < num_points; ++i)
    near_ipoints[i] = ptr_slist_new();
  while (iter < prob->max_iters)
  {
    // Assemble our points into a point set so that we can easily 
    // perform nearest-neighbor searches.
    kd_tree_t* tree = kd_tree_new(points, num_points);

    // Choose q points from within the domain according to the density 
    // function, and organize them into Voronoi regions of the points
    // in our point set.
    generate_random_points(prob->rng, density, bounding_box, prob->num_samples, samples);
    for (int j = 0; j < prob->num_samples; ++j)
    {
      int i = kd_tree_nearest(tree, &samples[j]);
      ptr_slist_append(near_ipoints[i], &samples[j]);
    }

    // Now we correct the generator positions.
    for (int i = 0; i < num_points; ++i)
    {
      ptr_slist_t* nearest = near_ipoints[i];
      correct_generator(alpha1, beta1, alpha2, beta2, nearest, &points[i], &ji[i]);
      ptr_slist_clear(nearest);

      // If we are outside the boundary (or within the minimum distance), 
      // project to the boundary.
      if (boundary != NULL)
      {
        real_t D;
        sp_func_eval(boundary, &points[i], &D);
        if ((fabs(D) < prob->min_dist) || (D > 0.0))
        {
          project_point_to_boundary(boundary, &points[i]);
          sp_func_eval(boundary, &points[i], &D);
          ASSERT(fabs(D) < 1e-12);
        }
      }
    }

    ++iter;
    kd_tree_free(tree);
  }

  // Clean up.
  for (int i = 0; i < num_points; ++i)
    ptr_slist_free(near_ipoints[i]);
}

cvt_gen_dist_t* prob_cvt_gen_dist_new(int (*random_gen)(), 
                                      int num_samples, 
                                      real_t alpha, 
                                      real_t beta, 
                                      real_t min_dist,
                                      int max_iters)
{
  ASSERT(num_samples > 0);
  ASSERT(alpha >= 0.0);
  ASSERT(alpha <= 1.0);
  ASSERT(beta >= 0.0);
  ASSERT(beta <= 1.0);
  ASSERT(min_dist >= 0.0);
  ASSERT(max_iters > 0);
  prob_cvt_gen_dist_t* prob = malloc(sizeof(prob_cvt_gen_dist_t));
  prob->num_samples = num_samples;
  prob->alpha1 = alpha;
  prob->alpha2 = 1.0 - alpha;
  prob->beta1 = beta;
  prob->beta2 = 1.0 - beta;
  prob->rng = random_gen;
  prob->min_dist = min_dist;
  prob->max_iters = max_iters;
  cvt_gen_dist_vtable vtable = {.iterate = prob_cvt_gen_dist_iterate, .dtor = free};
  return cvt_gen_dist_new("Probabilistic CVT generator distribution", prob, vtable);
}


