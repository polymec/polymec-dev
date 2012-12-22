#include <gc/gc.h>
#include "core/slist.h"
#include "core/point_set.h"
#include "geometry/prob_cvt_gen.h"

#ifdef __cplusplus
extern "C" {
#endif

struct prob_cvt_gen_t
{
  int q; // Number of sample points.
  double alpha1, beta1, alpha2, beta2; // Algorithm coefficients.
  long (*rng)(); // Random number generator.
};

prob_cvt_gen_t* prob_cvt_gen_new(long (*random_gen)(), int num_samples, double alpha, double beta)
{
  ASSERT(num_samples > 0);
  ASSERT(alpha >= 0.0);
  ASSERT(alpha <= 1.0);
  ASSERT(beta >= 0.0);
  ASSERT(beta <= 1.0);
  prob_cvt_gen_t* prob = GC_MALLOC(sizeof(prob_cvt_gen_t));
  prob->q = num_samples;
  prob->alpha1 = alpha;
  prob->alpha2 = 1.0 - alpha;
  prob->beta1 = beta;
  prob->beta2 = 1.0 - beta;
  prob->rng = random_gen;
  return prob;
}

// Termination critierion class.
struct prob_cvt_gen_term_t
{
  // A description of the termination condition.
  char* description;
  // Context pointer.
  void* context;
  // Termination critierion. Returns true if the iteration is to be terminated.
  bool (*terminate)(void*, point_t*, int, int);
  // Context destructor.
  void (*dtor)(void*);
};

// Destructor for termination critierion.
static void prob_cvt_gen_term_free(void* ctx, void* dummy)
{
  prob_cvt_gen_term_t* prob = (prob_cvt_gen_term_t*)ctx;
  free(prob->description);
  if ((prob->context != NULL) && (prob->dtor != NULL))
    prob->dtor(prob->context);
}

// Generic constructor for termination critierion.
static prob_cvt_gen_term_t* prob_cvt_gen_term_new(const char* description,
                                                  void* context,
                                                  bool (*terminate)(void*, point_t*, int, int),
                                                  void (*dtor)(void*))
{
  ASSERT(terminate != NULL);

  prob_cvt_gen_term_t* term = GC_MALLOC(sizeof(prob_cvt_gen_term_t));
  term->description = strdup(description);
  term->context = context;
  term->terminate = terminate;
  term->dtor = dtor;
  GC_register_finalizer(term, &prob_cvt_gen_term_free, term, NULL, NULL);
  return term;
}

static bool iteration_terminated(prob_cvt_gen_term_t* term, point_t* points, int num_points, int iteration)
{
  return term->terminate(term->context, points, num_points, iteration);
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

static void iterate(prob_cvt_gen_t* prob, 
                    sp_func_t* density,
                    sp_func_t* boundary,
                    bbox_t* bounding_box,
                    prob_cvt_gen_term_t* termination,
                    point_t* points, 
                    int num_points,
                    bool project_to_boundary)
{
  ASSERT(density != NULL);
  ASSERT(bounding_box != NULL);
  ASSERT(termination != NULL);
  ASSERT(points != NULL);
  ASSERT(num_points > 0);

  int iter = 0;
  int j[num_points];

  // Set ji to 1 for all i.
  for (int i = 0; i < num_points; ++i)
    j[i] = 1;

  double alpha1 = prob->alpha1;
  double beta1 = prob->beta1;
  double alpha2 = prob->alpha2;
  double beta2 = prob->beta2;

  // Iterate until termination.
  point_set_t* pset = point_set_new();
  while (!iteration_terminated(termination, points, num_points, iter))
  {
    // Assemble our points into a point set so that we can easily 
    // perform nearest-neighbor searches.
    point_set_clear(pset);
    for (int i = 0; i < num_points; ++i)
      point_set_insert(pset, &points[i], i);

    // Choose q points from within the domain according to the density 
    // function.
    point_t samples[prob->q];
    choose_sample_points(prob->rng, density, boundary, bounding_box, samples, prob->q);

    // Now organize the sample points into Voronoi regions of the points
    // in our point set.
    ptr_slist_t* near_points[num_points];
    for (int i = 0; i < num_points; ++i)
      near_points[i] = ptr_slist_new();
    for (int j = 0; j < prob->q; ++j)
    {
      int i = point_set_nearest(pset, &samples[j]);
      ptr_slist_append(near_points[i], &samples[j]);
    }

    // Now we correct the generator positions.
    for (int i = 0; i < num_points; ++i)
    {
      ptr_slist_t* nearest = near_points[i];

      point_t* zi = &points[i];
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
        int ji = j[i];
        zi->x = ((alpha1*ji + beta1)*zi->x + (alpha2*ji + beta2)*ui.x) / (ji + 1.0);
        zi->y = ((alpha1*ji + beta1)*zi->y + (alpha2*ji + beta2)*ui.y) / (ji + 1.0);
        zi->z = ((alpha1*ji + beta1)*zi->z + (alpha2*ji + beta2)*ui.z) / (ji + 1.0);

        // Project to the boundary if necessary.
        if (project_to_boundary)
        {
          ASSERT(sp_func_has_deriv(boundary, 1));
          double D, grad_D[3];
          sp_func_eval(boundary, zi, &D);
          sp_func_eval_deriv(boundary, 1, zi, grad_D);
          zi->x -= D * grad_D[0];
          zi->y -= D * grad_D[1];
          zi->z -= D * grad_D[2];
        }

        // Increment ji.
        ++(j[i]);
      }

      // Clean up.
      ptr_slist_free(nearest);
    }

    ++iter;
  }

  point_set_free(pset);
}

void prob_cvt_gen_iterate(prob_cvt_gen_t* prob, 
                          sp_func_t* density,
                          sp_func_t* boundary,
                          bbox_t* bounding_box,
                          prob_cvt_gen_term_t* termination,
                          point_t* points, 
                          int num_points)
{
  iterate(prob, density, boundary, bounding_box, termination, 
          points, num_points, false);
}

void prob_cvt_gen_iterate_on_boundary(prob_cvt_gen_t* prob, 
                                      sp_func_t* density,
                                      sp_func_t* boundary,
                                      bbox_t* bounding_box,
                                      prob_cvt_gen_term_t* termination,
                                      point_t* boundary_points, 
                                      int num_boundary_points)
{
  iterate(prob, density, boundary, bounding_box, termination, 
          boundary_points, num_boundary_points, true);
}

// Termination after max_iter steps.
typedef struct 
{
  int max_iter;
} max_iter_term_t;

static bool max_iter_terminate(void* context, point_t* points, int num_points, int iter)
{
  max_iter_term_t* m = (max_iter_term_t*)context;
  return (iter >= m->max_iter);
}

prob_cvt_gen_term_t* terminate_prob_cvt_at_iter(int max_iter)
{
  ASSERT(max_iter > 0);
  max_iter_term_t* context = malloc(sizeof(max_iter_term_t));
  context->max_iter = max_iter;
  char descr[1024];
  sprintf(descr, "Termination after %d iteration(s)", max_iter);
  return prob_cvt_gen_term_new((const char*)descr, context, max_iter_terminate, free);
}

#ifdef __cplusplus
}
#endif

