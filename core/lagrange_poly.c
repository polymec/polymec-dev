#include <gc/gc.h>
#include "core/lagrange_poly.h"

#ifdef __cplusplus
extern "C" {
#endif

struct lagrange_poly_t 
{
  int order;
  double *points, *basis, *weights;
};

static void lagrange_poly_free(void* ctx, void* dummy)
{
  lagrange_poly_t* poly = ctx;
  free(poly->weights);
  free(poly->basis);
  free(poly->points);
}

lagrange_poly_t* lagrange_poly_new(int order)
{
  ASSERT(order >= 1);
  lagrange_poly_t* poly = GC_MALLOC(sizeof(lagrange_poly_t));
  poly->order = order;
  poly->points = malloc(sizeof(double)*(order+1));
  poly->basis = malloc(sizeof(double)*(order+1));
  poly->weights = malloc(sizeof(double)*(order+1));
  GC_register_finalizer(poly, lagrange_poly_free, poly, NULL, NULL);
  return poly;
}

int lagrange_poly_order(lagrange_poly_t* poly)
{
  return poly->order;
}

void lagrange_poly_set_points(lagrange_poly_t* poly, double* points)
{
  // Set the interpolation points.
  memcpy(poly->points, points, sizeof(double)*(poly->order + 1));

  // Set the barycentric weights for quick evaluation.
  for (int j = 0; j < poly->order + 1; ++j)
  {
    double denom = 1.0;
    for (int i = 0; i < poly->order + 1; ++i)
    {
      if (i != j)
        denom *= (poly->points[j] - poly->points[i]);
    }
    poly->weights[j] = 1.0 / denom;
  }
}

void lagrange_poly_evaluate_basis(lagrange_poly_t* poly, double x, double* basis)
{
  for (int m = 0; m < poly->order+1; ++m)
  {
    basis[m] = 1.0;
    for (int j = 0; j < poly->order+1; ++j)
    {
      if (j != m)
        basis[m] *= (x - poly->points[m]) / (poly->points[j] - poly->points[m]);
    }
  }
}

double lagrange_poly_value(lagrange_poly_t* poly, double x, double* values)
{
  // Evaluate the polynomial using the barycentric representation.
  double lx = 1.0, sum = 0.0;
  for (int i = 0; i < poly->order+1; ++i)
  {
    double term = x - poly->points[i];
    if (term == 0.0)
      return values[i];
    lx *= term;
    sum += poly->weights[i] * values[i] / term;
  }
  return lx * sum;
}

void lagrange_poly_evaluate_basis_deriv(lagrange_poly_t* poly, int p, double x, double* basis_deriv)
{
  ASSERT(p >= 0);
  ASSERT(p <= 2);

  if (p > poly->order) 
  {
    memset(basis_deriv, 0, sizeof(double)*(poly->order + 1));
    return;
  }

  if (p == 0)
    lagrange_poly_evaluate_basis(poly, x, basis_deriv);
  else if (p == 1)
  {
    for (int m = 0; m < poly->order+1; ++m)
    {
      basis_deriv[m] = 0.0;
      for (int i = 0; i < poly->order+1; ++i)
      {
        if (i != m)
        {
          double prod = 1.0;
          for (int j = 0; j < poly->order+1; ++j)
          {
            if ((j != m) && (j != i))
              prod *= (x - poly->points[m]) / (poly->points[j] - poly->points[m]);
          }
          basis_deriv[m] += 1.0 / (poly->points[m] - poly->points[i]) * prod;
        }
      }
    }
  }
  else
  {
    ASSERT(p == 2);
    for (int m = 0; m < poly->order+1; ++m)
    {
      basis_deriv[m] = 0.0;
      for (int i = 0; i < poly->order+1; ++i)
      {
        if (i != m)
        {
          double sum = 0.0;
          for (int n = 0; n < poly->order+1; ++n)
          {
            double prod = 1.0;
            if (n != m)
            {
              for (int j = 0; j < poly->order+1; ++j)
              {
                if ((j != m) && (j != i) && (j != n))
                  prod *= (x - poly->points[j]) / (poly->points[m] - poly->points[j]);
              }
              sum += 1.0 / (poly->points[m] - poly->points[n]) * prod;
            }
          }
          basis_deriv[m] += 1.0 / (poly->points[m] - poly->points[i]) * sum;
        }
      }
    }
  }
}

double lagrange_poly_deriv(lagrange_poly_t* poly, int p, double x, double* values)
{
  // Evaluate the basis derivative at those points.
  lagrange_poly_evaluate_basis_deriv(poly, p, x, poly->basis);

  // Evaluate the polynomial.
  double val = 0.0;
  for (int j = 0; j < poly->order+1; ++j)
    val += values[j] * poly->basis[j];

  return val;
}

struct tensor_lagrange_poly_t 
{
  lagrange_poly_t* polys[3];
};

static void tensor_lagrange_poly_free(void* ctx, void* dummy)
{
  tensor_lagrange_poly_t* poly = ctx;
  for (int d = 0; d < 3; ++d)
    poly->polys[d] = NULL;
}

tensor_lagrange_poly_t* tensor_lagrange_poly_new(int order)
{
  tensor_lagrange_poly_t* poly = GC_MALLOC(sizeof(tensor_lagrange_poly_t));
  for (int d = 0; d < 3; ++d)
    poly->polys[d] = lagrange_poly_new(order);
  GC_register_finalizer(poly, tensor_lagrange_poly_free, poly, NULL, NULL);
  return poly;
}

int tensor_lagrange_poly_order(tensor_lagrange_poly_t* poly)
{
  return poly->polys[0]->order;
}

void tensor_lagrange_poly_set_points(tensor_lagrange_poly_t* poly, double* xs, double* ys, double* zs)
{
  lagrange_poly_set_points(poly->polys[0], xs);
  lagrange_poly_set_points(poly->polys[1], ys);
  lagrange_poly_set_points(poly->polys[2], zs);
}

double tensor_lagrange_poly_value(tensor_lagrange_poly_t* poly, point_t* x, double* values)
{
  double L1 = lagrange_poly_value(poly->polys[0], x->x, values);
  double L2 = lagrange_poly_value(poly->polys[1], x->y, values);
  double L3 = lagrange_poly_value(poly->polys[2], x->z, values);
  return L1*L2*L3;
}

double tensor_lagrange_poly_deriv(tensor_lagrange_poly_t* poly, int p, point_t* x, double* values)
{
  ASSERT(p >= 0);
  ASSERT(p <= 2);
  double L1 = lagrange_poly_value(poly->polys[0], x->x, values);
  double L2 = lagrange_poly_value(poly->polys[1], x->y, values);
  double L3 = lagrange_poly_value(poly->polys[2], x->z, values);
  if (p == 0)
    return L1*L2*L3;
  else if (p == 1)
  {
    double dL1 = lagrange_poly_deriv(poly->polys[0], 1, x->x, values);
    double dL2 = lagrange_poly_deriv(poly->polys[1], 1, x->y, values);
    double dL3 = lagrange_poly_deriv(poly->polys[2], 1, x->z, values);
    return dL1*L2*L3 + L1*dL2*L3 + L1*L2*dL3;
  }
  else
  {
    ASSERT(p == 2);
    double dL1 = lagrange_poly_deriv(poly->polys[0], 1, x->x, values);
    double dL2 = lagrange_poly_deriv(poly->polys[1], 1, x->y, values);
    double dL3 = lagrange_poly_deriv(poly->polys[2], 1, x->z, values);
    double d2L1 = lagrange_poly_deriv(poly->polys[0], 2, x->x, values);
    double d2L2 = lagrange_poly_deriv(poly->polys[1], 2, x->y, values);
    double d2L3 = lagrange_poly_deriv(poly->polys[2], 2, x->z, values);
    return d2L1*L2*L3 + dL1*dL2*L3 + dL1*L2*dL3 + 
           L1*d2L2*L3 + dL1*dL2*L3 + L1*dL2*dL3 + 
           L1*L2*d2L3 + dL1*L2*d2L3 + L1*dL2*dL3;
  }
}

#ifdef __cplusplus
}
#endif

