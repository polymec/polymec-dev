#include <gc/gc.h>
#include "geometry/polygon.h"
#include "core/slist.h"

#ifdef __cplusplus
extern "C" {
#endif

struct polygon_t 
{
  point_t* vertices;
  int num_vertices;
  double area;
  vector_t normal;
};

static void polygon_free(void* ctx, void* dummy)
{
  polygon_t* poly = ctx;
  free(poly->vertices);
}

static void polygon_compute_area(polygon_t* poly)
{
  // Compute the area using the fan algorithm.
  poly->area = 0.0;
  vector_t A, B;
  for (int j = 1; j < poly->num_vertices - 1; ++j)
  {
    // Form a triangle from vertex 0, vertex j, and vertex j+1.
    point_displacement(&poly->vertices[0], &poly->vertices[j], &A);
    point_displacement(&poly->vertices[0], &poly->vertices[j+1], &B);
    poly->area += vector_cross_mag(&A, &B);
  }
}

static void polygon_compute_normal(polygon_t* poly)
{
  vector_t A, B;
  point_displacement(&poly->vertices[0], &poly->vertices[1], &A);
  point_displacement(&poly->vertices[0], &poly->vertices[2], &B);
  vector_cross(&A, &B, &poly->normal);
  ASSERT(vector_mag(&poly->normal) > 0.0);
  vector_normalize(&poly->normal);
}

polygon_t* polygon_new(point_t* vertices, int num_vertices)
{
  ASSERT(vertices != NULL);
  ASSERT(num_vertices >= 3);
  polygon_t* poly = GC_MALLOC(sizeof(polygon_t));
  poly->vertices = malloc(sizeof(point_t)*num_vertices);
  memcpy(poly->vertices, vertices, sizeof(point_t)*num_vertices);
  poly->num_vertices = num_vertices;
  polygon_compute_area(poly);
  polygon_compute_normal(poly);
  GC_register_finalizer(poly, polygon_free, poly, NULL, NULL);
  return poly;
}

polygon_t* polygon_giftwrap(point_t* points, int num_points)
{
  ASSERT(num_points > 2);

  int count = 0;
  // FIXME

  point_t vertices[count];
  memcpy(vertices, points, sizeof(point_t)*count);
  return polygon_new(points, count);
}

int polygon_num_vertices(polygon_t* poly)
{
  return poly->num_vertices;
}

bool polygon_next_vertex(polygon_t* poly, int* pos, point_t** vertex)
{
  if (*pos >= poly->num_vertices) 
    return false;
  *vertex = &poly->vertices[*pos];
  ++(*pos);
  return true;
}

double polygon_area(polygon_t* poly)
{
  return poly->area;
}

polygon_t* polygon_clone(polygon_t* poly)
{
  return polygon_new(poly->vertices, poly->num_vertices);
}

void polygon_clip(polygon_t* poly, polygon_t* other)
{
  // FIXME
}

#ifdef __cplusplus
}
#endif

