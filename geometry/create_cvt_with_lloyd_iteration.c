#include "geometry/create_cvt_with_lloyd_iteration.h"

// Here we define a Lloyd CVT iterator and its functionality.
typedef struct
{
  point_t* centroids;
  int num_iterations, num_mobile_points;
} lloyd_cvt_iter_t;

static void lloyd_move_points(void* context, point_t* mobile_generators, int num_mobile_generators)
{
  lloyd_cvt_iter_t* lloyd = context;

  // Move each of the mobile generators to their centroids.
  memcpy(mobile_generators, lloyd->centroids, sizeof(point_t) * num_mobile_generators);
}

static bool lloyd_is_finished(void* context, mesh_t* mesh, int iteration)
{
  lloyd_cvt_iter_t* lloyd = context;

  // Jot down the centroids of the cells in the mesh.
  if (lloyd->centroids == NULL)
  for (int c = lloyd->num_mobile_points; c < mesh->num_cells; ++c)
    lloyd->centroids[c-lloyd->num_mobile_points] = mesh->cells[c].center;

  // We're finished if we reach the desired number of iterations.
  return (iteration >= lloyd->num_iterations); 
}

static void lloyd_dtor(void* context)
{
  lloyd_cvt_iter_t* lloyd = context;
  if (lloyd->centroids != NULL)
    free(lloyd->centroids);
  free(lloyd);
}

mesh_t* create_cvt_with_lloyd_iteration(point_t* stationary_generators, int num_stationary_generators, 
                                        point_t* mobile_generators, int num_mobile_generators,
                                        char** tag_names, int_array_t** tags, int num_tags,
                                        int num_iterations)
{
  lloyd_cvt_iter_t* lloyd = malloc(sizeof(lloyd_cvt_iter_t));
  lloyd->num_iterations = num_iterations;
  lloyd->num_mobile_points = num_mobile_generators;
  lloyd->centroids = malloc(sizeof(point_t) * num_mobile_generators);
  cvt_iterator_vtable vtable = {.move_points = lloyd_move_points,
                                .is_finished = lloyd_is_finished,
                                .dtor = lloyd_dtor};
  cvt_iterator_t* cvt_iter = cvt_iterator_new("Lloyd", lloyd, vtable);
  return create_cvt(stationary_generators, num_stationary_generators,
                    mobile_generators, num_mobile_generators,
                    tag_names, tags, num_tags, cvt_iter);
}

