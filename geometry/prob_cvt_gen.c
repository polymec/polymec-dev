#include <gc/gc.h>
#include "core/slist.h"
#include "geometry/prob_cvt_gen.h"

#ifdef __cplusplus
extern "C" {
#endif

struct prob_cvt_gen_t
{
  int q; // Number of sample points.
  double alpha1, beta1, alpha2, beta2; // Algorithm coefficients.
};

prob_cvt_gen_t* prob_cvt_gen_new(int num_samples, double alpha, double beta)
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

static void choose_sample_points(sp_func_t* density,
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
        p->x = (random()/RAND_MAX) * (bounding_box->x2 - bounding_box->x1) + bounding_box->x1;
        p->y = (random()/RAND_MAX) * (bounding_box->y2 - bounding_box->y1) + bounding_box->y1;
        p->z = (random()/RAND_MAX) * (bounding_box->z2 - bounding_box->z1) + bounding_box->z1;
        sp_func_eval(boundary, p, &Fp);
      }
      while (Fp >= 0.0);
    }
  }

  // Otherwise, just generate random points within the bounding box.
  else
  {
    for (int i = 0; i < num_points; ++i)
    {
      point_t* p = &points[i];
      p->x = (random()/RAND_MAX) * (bounding_box->x2 - bounding_box->x1) + bounding_box->x1;
      p->y = (random()/RAND_MAX) * (bounding_box->y2 - bounding_box->y1) + bounding_box->y1;
      p->z = (random()/RAND_MAX) * (bounding_box->z2 - bounding_box->z1) + bounding_box->z1;
    }
  }
}

static void find_near_points(point_t* z, point_t* samples, int num_samples, ptr_slist_t* near_points) 
{
  ASSERT(z != NULL);
  ASSERT(samples != NULL);
  ASSERT(num_samples > 0);
  ASSERT(near_points != NULL);
}

void prob_cvt_gen_iterate(prob_cvt_gen_t* prob, 
                          sp_func_t* density,
                          sp_func_t* boundary,
                          bbox_t* bounding_box,
                          prob_cvt_gen_term_t* termination,
                          point_t* points, 
                          int num_points)
{
  int iter = 0;
  int j[num_points];

  // Set ji to 1 for all i.
  for (int i = 0; i < num_points; ++i)
    j[i] = 1;

  double alpha1 = prob->alpha1;
  double beta1 = prob->beta1;
  double alpha2 = prob->alpha2;
  double beta2 = prob->beta2;

  // Set the random seed.
  // FIXME: This should probably be made more adjustable.
  unsigned int seed = 10;
  srandom(seed);

  // Iterate until termination.
  while (!iteration_terminated(termination, points, num_points, iter))
  {
    // Choose q points from within the domain according to the density 
    // function.
    point_t samples[prob->q];
    choose_sample_points(density, boundary, bounding_box, samples, prob->q);

    for (int i = 0; i < num_points; ++i)
    {
      // Find the sample points that are within the Voronoi region of 
      // the ith generator.
      ptr_slist_t* near_points = ptr_slist_new();
      point_t* zi = &points[i];
      find_near_points(zi, samples, prob->q, near_points);

      // Compute the average, ui, of the sample points in the Voronoi region.
      if (near_points->size > 0)
      {
        point_t ui = {.x = 0.0, .y = 0.0, .z = 0.0};
        ptr_slist_node_t* node = near_points->front;
        while (node != near_points->back)
        {
          point_t* pos = node->value;
          ui.x += pos->x;
          ui.y += pos->y;
          ui.z += pos->z;
        }
        ui.x /= near_points->size;
        ui.y /= near_points->size;
        ui.z /= near_points->size;

        // Correct the position of the generator.
        int ji = j[i];
        zi->x = ((alpha1*ji + beta1)*zi->x + (alpha2*ji + beta2)*ui.x) / (ji + 1.0);
        zi->y = ((alpha1*ji + beta1)*zi->y + (alpha2*ji + beta2)*ui.y) / (ji + 1.0);
        zi->z = ((alpha1*ji + beta1)*zi->z + (alpha2*ji + beta2)*ui.z) / (ji + 1.0);

        // Increment ji.
        ++(j[i]);
      }
    }

    ++iter;
  }
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

