#include <gc/gc.h>
#include "geometry/polygon.h"
#include "geometry/polygon2.h"
#include "geometry/plane.h"
#include "core/slist.h"

#ifdef __cplusplus
extern "C" {
#endif

struct polygon_t 
{
  point_t* vertices;
  int num_vertices;
  double area;
  point_t x0;
  vector_t normal;
  sp_func_t* plane;
};

static void polygon_free(void* ctx, void* dummy)
{
  polygon_t* poly = ctx;
  free(poly->vertices);
  poly->plane = NULL;
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
    poly->area += 0.5 * vector_cross_mag(&A, &B);
  }
}

static void compute_normal(point_t* points, int num_points, vector_t* n)
{
  vector_t A, B;
  point_displacement(&points[0], &points[1], &A);
  point_displacement(&points[0], &points[2], &B);
  vector_cross(&A, &B, n);
  ASSERT(vector_mag(n) > 0.0);
  vector_normalize(n);
}

static void compute_centroid(point_t* points, int num_points, point_t* x0)
{
  ASSERT(num_points > 0);
  x0->x = x0->y = x0->z = 0.0;
  for (int i = 0; i < num_points; ++i)
  {
    x0->x += points[i].x;
    x0->y += points[i].y;
    x0->z += points[i].z;
  }
  x0->x /= num_points;
  x0->y /= num_points;
  x0->z /= num_points;
}

static void polygon_compute_plane(polygon_t* poly)
{
  compute_normal(poly->vertices, poly->num_vertices, &poly->normal);
  compute_centroid(poly->vertices, poly->num_vertices, &poly->x0);
  poly->plane = plane_new(&poly->normal, &poly->x0);
}

polygon_t* polygon_new(point_t* vertices, int num_vertices)
{
  ASSERT(vertices != NULL);
  ASSERT(num_vertices >= 3);
  // FIXME: Check that all points are coplanar.
  polygon_t* poly = GC_MALLOC(sizeof(polygon_t));
  poly->vertices = malloc(sizeof(point_t)*num_vertices);
  memcpy(poly->vertices, vertices, sizeof(point_t)*num_vertices);
  poly->num_vertices = num_vertices;
  polygon_compute_area(poly);
  polygon_compute_plane(poly);
  GC_register_finalizer(poly, polygon_free, poly, NULL, NULL);
  return poly;
}

polygon_t* polygon_giftwrap(point_t* points, int num_points)
{
  ASSERT(num_points > 2);
  // FIXME: Check that all points are coplanar.

  // Find the plane of the polygon.
  vector_t normal;
  compute_normal(points, num_points, &normal);
  point_t x0;
  compute_centroid(points, num_points, &x0);
  sp_func_t* plane = plane_new(&normal, &x0);

  // Do the gift-wrapping in 2D.
  point2_t pts[num_points];
  for (int i = 0; i < num_points; ++i)
    plane_project(plane, &points[i], &pts[i]);
  polygon2_t* poly2 = polygon2_giftwrap(pts, num_points);

  // Re-embed the resulting vertices in 3D.
  int num_vertices = polygon2_num_vertices(poly2);
  point_t vertices[num_vertices];
  int pos = 0, offset = 0;
  point2_t* vtx;
  while (polygon2_next_vertex(poly2, &pos, &vtx))
    plane_embed(plane, vtx, &vertices[offset++]);

  // Clean up.
  poly2 = NULL;
  plane = NULL;

  return polygon_new(vertices, num_vertices);
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
  // Do the clipping in 2D.
  point2_t pts[poly->num_vertices];
  for (int i = 0; i < poly->num_vertices; ++i)
    plane_project(poly->plane, &poly->vertices[i], &pts[i]);
  polygon2_t* poly2 = polygon2_new(pts, poly->num_vertices);

  point2_t other_pts[other->num_vertices];
  for (int i = 0; i < other->num_vertices; ++i)
    plane_project(other->plane, &other->vertices[i], &other_pts[i]);
  polygon2_t* other2 = polygon2_new(other_pts, other->num_vertices);

  polygon2_clip(poly2, other2);

  // Now re-embed the vertices.
  int num_vertices = polygon2_num_vertices(poly2);
  if (poly->num_vertices < num_vertices)
    poly->vertices = realloc(poly->vertices, sizeof(point_t)*num_vertices);
  poly->num_vertices = num_vertices;
  int pos = 0, offset = 0;
  point2_t* vtx;
  while (polygon2_next_vertex(poly2, &pos, &vtx))
    plane_embed(poly->plane, vtx, &poly->vertices[offset++]);

  // Recompute geometry (plane remains unchanged).
  compute_centroid(poly->vertices, poly->num_vertices, &poly->x0);

  // Clean up.
  poly2 = NULL;
  other2 = NULL;
}

#ifdef __cplusplus
}
#endif

