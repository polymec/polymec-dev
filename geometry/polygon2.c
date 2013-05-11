#include <gc/gc.h>
#include "geometry/polygon2.h"
#include "core/slist.h"

struct polygon2_t 
{
  point2_t* vertices;
  int num_vertices;
  int* ordering;
  double area;
};

static void polygon2_free(void* ctx, void* dummy)
{
  polygon2_t* poly = ctx;
  free(poly->vertices);
  free(poly->ordering);
}

static void polygon2_compute_area(polygon2_t* poly)
{
  // Compute the area using the fan algorithm.
  poly->area = 0.0;
  int I = poly->ordering[0];
  for (int j = 1; j < poly->num_vertices - 1; ++j)
  {
    // Form a triangle from vertex 0, vertex j, and vertex j+1.
    int J = poly->ordering[j];
    int K = poly->ordering[j+1];
    poly->area += triangle_area(&poly->vertices[I], &poly->vertices[J], &poly->vertices[K]);
  }
}

polygon2_t* polygon2_new(point2_t* vertices, int num_vertices)
{
  int ordering[num_vertices];
  for (int i = 0; i < num_vertices; ++i)
    ordering[i] = i;
  return polygon2_new_with_ordering(vertices, ordering, num_vertices);
}

polygon2_t* polygon2_new_with_ordering(point2_t* points, int* ordering, int num_points)
{
  ASSERT(points != NULL);
  ASSERT(num_points >= 3);
  polygon2_t* poly = GC_MALLOC(sizeof(polygon2_t));
  poly->vertices = malloc(sizeof(point2_t)*num_points);
  memcpy(poly->vertices, points, sizeof(point2_t)*num_points);
  poly->num_vertices = num_points;
  poly->ordering = malloc(sizeof(int)*num_points);
  memcpy(poly->ordering, ordering, sizeof(int)*num_points);
  polygon2_compute_area(poly);
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

  return polygon2_new_with_ordering(points, indices, count);
}

int polygon2_num_vertices(polygon2_t* poly)
{
  return poly->num_vertices;
}

int* polygon2_ordering(polygon2_t* poly)
{
  return poly->ordering;
}

bool polygon2_next_vertex(polygon2_t* poly, int* pos, point2_t** vertex)
{
  if (*pos >= poly->num_vertices) 
    return false;
  *vertex = &poly->vertices[poly->ordering[*pos]];
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

static char parallel_int(point2_t* a, point2_t* b, point2_t* c, 
                         point2_t* d, point2_t* p)
{
  if (!point2s_are_colinear(a, b, c))
    return '0';
  if (point2_is_between(c, a, b))
  {
    point2_copy(p, c);
    return 'e';
  }
  if (point2_is_between(d, a, b))
  {
    point2_copy(p, d);
    return 'e';
  }
  if (point2_is_between(a, c, d))
  {
    point2_copy(p, a);
    return 'e';
  }
  if (point2_is_between(b, c, d))
  {
    point2_copy(p, b);
    return 'e';
  }
  return '0';
}

static char seg_set_int(point2_t* a, point2_t* b, point2_t* c, 
                        point2_t* d, point2_t* p)
{
  char code;

  double denom = a->x * (d->y - c->y) + 
                 b->x * (c->y - d->y) + 
                 d->x * (b->y - a->y) + 
                 c->x * (a->y - b->y);

  // If denom is zero, segments are parallel.
  if (denom == 0.0)
    return parallel_int(a, b, c, d, p);

  double num = a->x * (d->y - c->y) + 
               c->x * (a->y - d->y) + 
               d->x * (c->y - a->y);
  if ((num == 0.0) || (num == denom)) 
    code = 'v';
  double s = num / denom;

  num = -(a->x * (c->y - b->y) + 
          b->x * (a->y - c->y) + 
          c->x * (b->y - a->y));
  if ((num == 0.0) || (num == denom))
    code = 'v';
  double t = num / denom;

  if ((0.0 < s) && (s < 1.0) &&
      (0.0 < t) && (t < 1.0))
    code = '1';
  else if ((0.0 > s) || (s > 1.0) || 
           (0.0 > t) || (t > 1.0))
    code = '0';

  p->x = a->x + s * (b->x - a->x);
  p->y = a->y + s * (b->y - a->y);

  return code;
}

typedef enum
{
  UNKNOWN,
  PIN,
  QIN
} polygon2_inout_t;

static polygon2_inout_t in_out(point2_t* p, polygon2_inout_t inflag, int aHB, int bHA, double_slist_t* xlist, double_slist_t* ylist)
{
  double_slist_append(xlist, p->x);
  double_slist_append(ylist, p->y);
  if (aHB > 0)
    return PIN;
  else if (bHA > 0)
    return QIN;
  else 
    return inflag;
}

static int advance(int a, int* aa, int n, bool inside, point2_t* v, double_slist_t* xlist, double_slist_t* ylist)
{
  if (inside)
  {
    double_slist_append(xlist, v->x);
    double_slist_append(ylist, v->y);
  }
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
  static point2_t origin = {.x = 0.0, .y = 0.0};
  point2_t p0; // First point.

  // We'll store the coordinates of the vertices of the 
  // clipped polygon here.
  double_slist_t* xlist = double_slist_new();
  double_slist_t* ylist = double_slist_new();

  do
  {
    // Computations of key variables.
    int a1 = (a + n - 1) % n;
    int b1 = (b + m - 1) % m;

    int aa = poly->ordering[a];
    int aa1 = poly->ordering[a1];
    int bb = poly->ordering[b];
    int bb1 = poly->ordering[b1];

    point2_t* Pa = &poly->vertices[aa];
    point2_t* Pa1 = &poly->vertices[aa1];
    point2_t* Qb = &other->vertices[bb];
    point2_t* Qb1 = &other->vertices[bb1];

    point2_t A, B;
    A.x = Pa1->x - Pa->x;
    A.y = Pa1->y - Pa->y;
    B.x = Qb1->x - Qb->x;
    B.y = Qb1->y - Qb->y;

    int cross = SIGN(triangle_area(&origin, &A, &B));
    int aHB   = SIGN(triangle_area(Qb1, Qb, Pa));
    int bHA   = SIGN(triangle_area(Pa1, Pa, Qb));

    // If A and B intersect, update inflag.
    point2_t p;
    char code = seg_set_int(Pa1, Pa, Qb1, Qb, &p);
    if ((code == '1') || (code == 'v')) 
    {
      if ((inflag == UNKNOWN) && (double_slist_empty(xlist)))
      {
        aa = ba = 0;
        p0.x = p.x;
        p0.y = p.y;
        double_slist_append(xlist, p0.x);
        double_slist_append(ylist, p0.y);
      }
      inflag = in_out(&p, inflag, aHB, bHA, xlist, ylist);
    }

    // Advance rules.
    else if (cross >= 0)
    {
      if (bHA > 0)
        a = advance(a, &aa, n, (inflag == PIN), Pa, xlist, ylist);
      else
        b = advance(b, &ba, m, (inflag == QIN), Qb, xlist, ylist);
    }
    else // if (cross < 0)
    {
      if (aHB > 0)
        b = advance(b, &ba, m, (inflag == QIN), Qb, xlist, ylist);
      else
        a = advance(a, &aa, n, (inflag == PIN), Pa, xlist, ylist);
    }
  }
  while (((aa < n) || (ba < m)) && (aa < 2*n) && (ba < 2*m));

  // Replace this polygon with its clipped version.
  ASSERT(xlist->size == ylist->size);
  if (xlist->size > poly->num_vertices)
    realloc(poly->vertices, sizeof(point2_t)*xlist->size);
  poly->num_vertices = xlist->size;
  for (int i = 0; i < xlist->size; ++i)
  {
    // FIXME: Verify!
    poly->vertices[i].x = double_slist_pop(xlist);
    poly->vertices[i].y = double_slist_pop(ylist);
  }
  polygon2_compute_area(poly);

  // Clean up.
  double_slist_free(xlist);
  double_slist_free(ylist);
}

