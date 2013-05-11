#include <gc/gc.h>
#include "core/lagrange_poly.h"

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
              prod *= (x - poly->points[j]) / (poly->points[m] - poly->points[j]);
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
            if ((n != i) && (n != m))
            {
              double prod = 1.0;
              for (int j = 0; j < poly->order+1; ++j)
              {
                if ((j != m) && (j != i) && (j != n))
                  prod *= (x - poly->points[j]) / (poly->points[m] - poly->points[j]);
              }
              sum += prod / (poly->points[m] - poly->points[n]);
            }
          }
          basis_deriv[m] += sum / (poly->points[m] - poly->points[i]);
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

