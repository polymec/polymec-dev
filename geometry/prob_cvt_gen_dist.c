#include "core/slist.h"
#include "core/point_set.h"
#include "geometry/prob_cvt_gen_dist.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct 
{
  int num_interior_samples, num_boundary_samples; // Number of sample points.
  double alpha1, beta1, alpha2, beta2; // Algorithm coefficients.
  long (*rng)(); // Random number generator.
  int max_iters; // Maximum number of iterations.
} prob_cvt_gen_dist_t;

static void project_point_to_boundary(sp_func_t* boundary, point_t* x)
{
  ASSERT(sp_func_has_deriv(boundary, 1));
  double D, grad_D[3];
  sp_func_eval(boundary, x, &D);
  sp_func_eval_deriv(boundary, 1, x, grad_D);
  vector_t normal = {.x = -grad_D[0], .y = -grad_D[1], .z = -grad_D[2]};
  vector_normalize(&normal);
//printf("%g %g %g (%g, %g, %g, %g) ->", x->x, x->y, x->z, D, normal.x, normal.y, normal.z);
  x->x += D * normal.x;
  x->y += D * normal.y;
  x->z += D * normal.z;
//printf("%g %g %g\n", x->x, x->y, x->z);
}

static void choose_sample_points(long (*rng)(),
                                 sp_func_t* density,
                                 sp_func_t* boundary,
                                 bbox_t* bounding_box,
                                 point_t* points,
                                 int num_points)
{
  ASSERT(density != NULL);
  ASSERT(bounding_box != NULL);
  ASSERT(points != NULL);
  ASSERT(num_points > 0);

  // Use a boundary function if given.
  if (boundary != NULL)
  {
    for (int i = 0; i < num_points; ++i)
    {
      point_t* p = &points[i];
      double Fp;
      do
      {
        point_randomize(p, rng, bounding_box);
        sp_func_eval(boundary, p, &Fp);
      }
      while (Fp >= 0.0);
    }
  }

  // Otherwise, just generate random points within the bounding box.
  else
  {
    for (int i = 0; i < num_points; ++i)
      point_randomize(&points[i], rng, bounding_box);
  }
}

static void correct_generator(double alpha1, double beta1, 
                              double alpha2, double beta2, 
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
                               point_t* interior_points, 
                               int num_interior_points,
                               point_t* boundary_points, 
                               int num_boundary_points)
{
  ASSERT(density != NULL);
  ASSERT(bounding_box != NULL);
  ASSERT(num_interior_points >= 0);
  ASSERT((interior_points != NULL) || (num_interior_points == 0));
  ASSERT(num_boundary_points >= 0);
  ASSERT((boundary_points != NULL) || (num_boundary_points == 0));

  prob_cvt_gen_dist_t* prob = context;

  int iter = 0;
  int ji[num_interior_points], jb[num_boundary_points];

  // Project the boundary points to the boundary.
  if ((num_boundary_points > 0) && (boundary != NULL))
  {
    for (int i = 0; i < num_boundary_points; ++i)
    {
      project_point_to_boundary(boundary, &boundary_points[i]);
      double D;
      sp_func_eval(boundary, &boundary_points[i], &D);
      if (fabs(D) > 1e-12)
        polymec_error("cvt_gen_dist_iterate: boundary projection yielded a non-zero distance (%g).\n Boundary is not a signed distance function!", D);
    }
  }

  // Set ji to 1 for all i.
  for (int i = 0; i < num_interior_points; ++i)
    ji[i] = 1;
  for (int i = 0; i < num_boundary_points; ++i)
    jb[i] = 1;

  double alpha1 = prob->alpha1;
  double beta1 = prob->beta1;
  double alpha2 = prob->alpha2;
  double beta2 = prob->beta2;

  // Iterate until termination.
  point_set_t* interior_pset = point_set_new();
  point_set_t* boundary_pset = point_set_new();
  point_t interior_samples[prob->num_interior_samples], 
          boundary_samples[prob->num_boundary_samples];
  ptr_slist_t* near_ipoints[num_interior_points];
  ptr_slist_t* near_bpoints[num_boundary_points];
  for (int i = 0; i < num_interior_points; ++i)
    near_ipoints[i] = ptr_slist_new();
  for (int i = 0; i < num_boundary_points; ++i)
    near_bpoints[i] = ptr_slist_new();
  while (iter < prob->max_iters)
  {
    // Assemble our points into a point set so that we can easily 
    // perform nearest-neighbor searches.
    point_set_clear(interior_pset);
    for (int i = 0; i < num_interior_points; ++i)
      point_set_insert(interior_pset, &interior_points[i], i);
    point_set_clear(boundary_pset);
    for (int i = 0; i < num_boundary_points; ++i)
      point_set_insert(boundary_pset, &boundary_points[i], i + num_interior_points);

    // Choose q points from within the domain according to the density 
    // function, and organize them into Voronoi regions of the points
    // in our point set.
    if (num_interior_points > 0)
    {
      choose_sample_points(prob->rng, density, boundary, bounding_box, interior_samples, prob->num_interior_samples);
      for (int j = 0; j < prob->num_interior_samples; ++j)
      {
        int i = point_set_nearest(interior_pset, &interior_samples[j]);
        ptr_slist_append(near_ipoints[i], &interior_samples[j]);
      }
    }
    if (num_boundary_points > 0)
    {
      choose_sample_points(prob->rng, density, boundary, bounding_box, boundary_samples, prob->num_boundary_samples);
      for (int i = 0; i < prob->num_boundary_samples; ++i)
      {
        project_point_to_boundary(boundary, &boundary_samples[i]);
        int j = point_set_nearest(boundary_pset, &boundary_samples[i]);
        ptr_slist_append(near_bpoints[j], &boundary_samples[i]);
      }
    }

    // Now we correct the generator positions.
    for (int i = 0; i < num_interior_points; ++i)
    {
      ptr_slist_t* nearest = near_ipoints[i];
      correct_generator(alpha1, beta1, alpha2, beta2, nearest, &interior_points[i], &ji[i]);
      ptr_slist_clear(nearest);
    }
    for (int i = 0; i < num_boundary_points; ++i)
    {
      ptr_slist_t* nearest = near_bpoints[i];
      correct_generator(alpha1, beta1, alpha2, beta2, nearest, &boundary_points[i], &jb[i]);

      // Make sure that the point is projected to the boundary.
      project_point_to_boundary(boundary, &boundary_points[i]);
      double D;
      sp_func_eval(boundary, &boundary_points[i], &D);
      ASSERT(fabs(D) < 1e-12);

      ptr_slist_clear(nearest);
    }

    ++iter;
  }

  // Clean up.
  for (int i = 0; i < num_interior_points; ++i)
    ptr_slist_free(near_ipoints[i]);
  for (int i = 0; i < num_boundary_points; ++i)
    ptr_slist_free(near_bpoints[i]);

  point_set_free(interior_pset);
  point_set_free(boundary_pset);
}

cvt_gen_dist_t* prob_cvt_gen_dist_new(long (*random_gen)(), 
                                      int num_interior_samples, 
                                      int num_boundary_samples,
                                      double alpha, 
                                      double beta, 
                                      int max_iters)
{
  ASSERT(num_interior_samples > 0);
  ASSERT(num_boundary_samples >= 0);
  ASSERT(alpha >= 0.0);
  ASSERT(alpha <= 1.0);
  ASSERT(beta >= 0.0);
  ASSERT(beta <= 1.0);
  ASSERT(max_iters > 0);
  prob_cvt_gen_dist_t* prob = malloc(sizeof(prob_cvt_gen_dist_t));
  prob->num_interior_samples = num_interior_samples;
  prob->num_boundary_samples = num_interior_samples;
  prob->alpha1 = alpha;
  prob->alpha2 = 1.0 - alpha;
  prob->beta1 = beta;
  prob->beta2 = 1.0 - beta;
  prob->rng = random_gen;
  prob->max_iters = max_iters;
  cvt_gen_dist_vtable vtable = {.iterate = prob_cvt_gen_dist_iterate, .dtor = free};
  return cvt_gen_dist_new("Probabilistic CVT generator distribution", prob, vtable);
}


#ifdef __cplusplus
}
#endif

