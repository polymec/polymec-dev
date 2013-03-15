#include <gc/gc.h>
#include "geometry/hexahedron.h"
#include "core/lagrange_poly.h"

#ifdef __cplusplus
extern "C" {
#endif

struct hexahedron_t 
{
  int order, num_points;
  lagrange_poly_t* poly;
};

static void hexahedron_free(void* ctx, void* dummy)
{
  hexahedron_t* hex = ctx;
  hex->poly = NULL;
}


hexahedron_t* hexahedron_new(int order)
{
  ASSERT(order >= 1);
  ASSERT(order <= 3);
  hexahedron_t* hex = GC_MALLOC(sizeof(hexahedron_t));
  hex->order = order;
  hex->num_points = pow(order+1, 3);
  hex->poly = lagrange_poly_new(order);

  // Set the points of the hex to those of the reference (logical) 
  // coordinate system.
  double points[order+1];
  double h = 1.0 / order;
  for (int i = 0; i < order+1; ++i)
    points[i] = i*h;
  lagrange_poly_set_points(hex->poly, points);

  GC_register_finalizer(hex, hexahedron_free, hex, NULL, NULL);

  return hex;
}

int hexahedron_order(hexahedron_t* hex)
{
  return hex->order;
}

int hexahedron_num_points(hexahedron_t* hex)
{
  return hex->num_points;
}

void hexahedron_get_points(hexahedron_t* hex, point_t* points)
{
  double h = 1.0 / hex->order;
  int offset = 0;
  for (int k = 0; k < hex->order+1; ++k)
  {
    for (int j = 0; j < hex->order+1; ++j)
    {
      for (int i = 0; i < hex->order+1; ++i, ++offset)
      {
        points[offset].z = k*h;
        points[offset].y = j*h;
        points[offset].z = i*h;
      }
    }
  }
}

void hexahedron_map(hexahedron_t* hex, point_t* points, point_t* xi, point_t* x)
{
  // Evaluate the Lagrange polynomials in each direction
  // and form their tensor product.
  int order = hex->order;
  double ones[order+1];
  for (int i = 0; i < order+1; ++i)
    ones[i] = 1.0;

  x->x = x->y = x->z = 0.0;
  int p = 0;
  for (int k = 0; k < order+1; ++k)
  {
    double Lk = lagrange_poly_value(hex->poly, xi->z, ones);
    for (int j = 0; j < order+1; ++j)
    {
      double Lj = lagrange_poly_value(hex->poly, xi->y, ones);
      for (int i = 0; i < order+1; ++i, ++p)
      {
        double Li = lagrange_poly_value(hex->poly, xi->x, ones);
        double LiLjLk = Li * Lj * Lk;
        x->x += points[p].x * LiLjLk;
        x->y += points[p].y * LiLjLk;
        x->z += points[p].z * LiLjLk;
      }
    }
  }
}

#ifdef __cplusplus
}
#endif

