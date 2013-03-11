#include <gc/gc.h>
#include "geometry/polygon2.h"

#ifdef __cplusplus
extern "C" {
#endif

struct polygon2_t 
{
  point2_t* vertices;
  int num_vertices;
  double area;
};

static void polygon2_free(void* ctx, void* dummy)
{
  polygon2_t* poly = ctx;
  free(poly->vertices);
}

static inline double triangle_area(point2_t* x1, point2_t* x2, point2_t* x3)
{
  return 0.5 * ((x2->x - x1->x) * (x3->y - x1->y) - (x3->x - x1->x) * (x2->y - x1->y));
}

polygon2_t* polygon2_new(point2_t* vertices, int num_vertices)
{
  ASSERT(vertices != NULL);
  ASSERT(num_vertices >= 3);
  polygon2_t* poly = GC_MALLOC(sizeof(polygon2_t));
  poly->vertices = malloc(sizeof(point2_t)*num_vertices);
  memcpy(poly->vertices, vertices, sizeof(point2_t)*num_vertices);
  poly->num_vertices = num_vertices;

  // Compute the area using the fan algorithm.
  poly->area = 0.0;
  for (int j = 1; j < num_vertices - 1; ++j)
  {
    // Form a triangle from vertex 0, vertex j, and vertex j+1.
    poly->area += triangle_area(&vertices[0], &vertices[j], &vertices[j+1]);
  }

  GC_register_finalizer(poly, polygon2_free, poly, NULL, NULL);
  return poly;
}

polygon2_t* polygon2_giftwrap(point2_t* points, int num_points)
{
  ASSERT(num_points > 2);

  int indices[num_points], count = 0;

  // Find the "lowest" point in the set.
  double ymin = FLT_MAX;
  int index0 = -1;
  for (int p = 0; p < num_points; ++p)
  {
    if (ymin > points[p].y)
    {
      ymin = points[p].y;
      index0 = p;
    }
  }

  // We start with this point and a horizontal angle.
  double theta_prev = 0.0;
  indices[count++] = index0;

  // Now start gift wrapping.
  int i = index0;
  do 
  {
    double dtheta_min = 2.0*M_PI;
    int j_min = -1;
    for (int j = 0; j < num_points; ++j)
    {
      if (j != i)
      {
        double dx = points[j].x - points[i].x,
               dy = points[j].y - points[i].y;
        double theta = atan2(dy, dx);
        double dtheta = theta - theta_prev;
        if (dtheta < 0.0)
          dtheta += 2.0*M_PI;
        if (dtheta_min > dtheta)
        {
          dtheta_min = dtheta;
          j_min = j;
        }
      }
    }
    if (j_min != index0)
      indices[count++] = j_min;
    theta_prev += dtheta_min;
    i = j_min;
  }
  while (i != index0);

  // The convex hull should be a polygon.
  ASSERT(count > 2);

  point2_t vertices[count];
  memcpy(vertices, points, sizeof(point2_t)*count);
  return polygon2_new(points, count);
}

int polygon2_num_vertices(polygon2_t* poly)
{
  return poly->num_vertices;
}

bool polygon2_next_vertex(polygon2_t* poly, int* pos, point2_t** vertex)
{
  if (*pos >= poly->num_vertices) 
    return false;
  *vertex = &poly->vertices[*pos];
  ++(*pos);
  return true;
}

double polygon2_area(polygon2_t* poly)
{
  return poly->area;
}

polygon2_t* polygon2_clone(polygon2_t* poly)
{
  return polygon2_new(poly->vertices, poly->num_vertices);
}

typedef enum
{
  UNKNOWN,
  IN,
  OUT
} polygon2_inout_t;

static int advance(int a, int* aa, int n, bool inside, point2_t* v)
{
  (*aa)++;
  return (a+1) % n;
}

void polygon2_clip(polygon2_t* poly, polygon2_t* other)
{
  // This implementation is based on that shown in Chapter 7.6.1 of 
  // Joseph O'Rourke's _Computational_Geometry_In_C_.

  int a = 0, b = 0, aa = 0, ba = 0;
  int n = poly->num_vertices, m = other->num_vertices;
  polygon2_inout_t inflag = UNKNOWN;
  bool first_point = true;
  static point2_t origin = {.x = 0.0, .y = 0.0};
  do
  {
    // Computations of key variables.
    int a1 = (a + n - 1) % n;
    int b1 = (b + m - 1) % m;

    point2_t A, B;
    A.x = poly->vertices[a1].x - poly->vertices[a].x;
    A.y = poly->vertices[a1].y - poly->vertices[a].y;
    B.x = poly->vertices[b1].x - poly->vertices[b].x;
    B.y = poly->vertices[b1].y - poly->vertices[b].y;

    int cross = SIGN(triangle_area(&origin, &A, &B));
    int aHB   = SIGN(triangle_area(&other->vertices[b1], &other->vertices[b], &poly->vertices[a]));
    int bHA   = SIGN(triangle_area(&poly->vertices[a1], &poly->vertices[a], &other->vertices[b]));

    // If A and B intersect, update inflag.
    // FIXME

    // Advance rules.
    // FIXME
  }
  while (((aa < n) || (ba < m)) && (aa < 2*n) && (ba < 2*m));
}

#ifdef __cplusplus
}
#endif

