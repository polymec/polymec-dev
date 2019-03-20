// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/lookup1.h"
#include "core/linear_algebra.h"

struct lookup1_t
{
  real_t x_min, x_max, dx;
  int num_values;
  real_t *values, *fit_data;
  lookup1_interp_t interp;
  real_t (*compute)(lookup1_t*, real_t);
};

static real_t lookup1_linear(lookup1_t* table, real_t x)
{
  // Identify the index of the y value to the left of x.
  real_t dx = table->dx;
  int index = (int)((x - table->x_min) / dx);
  ASSERT(index >= 0);
  ASSERT(index < table->num_values);
  real_t delta_x = x - index * dx;
  return table->values[index] + table->fit_data[index] * delta_x;
}

static real_t lookup1_quadratic(lookup1_t* table, real_t x)
{
  // Identify the index of the y value to the left of x.
  real_t dx = table->dx;
  int index = (int)((x - table->x_min) / dx);
  ASSERT(index >= 0);
  ASSERT(index < table->num_values);
  real_t delta_x = x - index * dx;
  return table->values[index] +
         table->fit_data[2*index] * delta_x +
         table->fit_data[2*index+1] * delta_x * delta_x;
}

lookup1_t* lookup1_new(real_t x_min, real_t x_max, int num_values,
                       real_t* y_values, lookup1_interp_t interpolation)
{
  ASSERT(x_min < x_max);
  ASSERT(num_values > 1);

  lookup1_t* table = polymec_malloc(sizeof(lookup1_t));
  table->x_min = x_min;
  table->x_max = x_max;
  table->num_values = num_values;
  table->dx = (x_max - x_min) / (num_values - 1);
  table->values = polymec_malloc(sizeof(real_t) * num_values);
  memcpy(table->values, y_values, sizeof(real_t) * num_values);
  table->interp = interpolation;
  if (table->interp == LOOKUP1_LINEAR)
  {
    // Use the linear lookup function.
    table->compute = lookup1_linear;

    // Fit data is dy/dx between every two points.
    table->fit_data = polymec_malloc(sizeof(real_t) * num_values);
    for (int i = 0; i < num_values-1; ++i)
    {
      real_t y1 = table->values[i];
      real_t y2 = table->values[i+1];
      table->fit_data[i] = (y2 - y1) / table->dx;
    }
    table->fit_data[num_values-1] = 0.0; // padding;
  }
  else
  {
    ASSERT(table->interp == LOOKUP1_QUADRATIC);
    table->compute = lookup1_quadratic;

    // Fit data is dy/dx and d2y/dx2 between every two points.
    table->fit_data = polymec_malloc(sizeof(real_t) * 2 * num_values);

    real_t dx = table->dx;

    // The first interval is a quadratic fit about the first of 3 points.
    {
      real_t y1 = table->values[0];
      real_t y2 = table->values[1];
      real_t y3 = table->values[2];

      real_t A[4] = {dx, 2*dx, dx*dx, 4*dx*dx};
      real_t B[2] = {y2-y1, y3-y1};
      real_t X[2];
      solve_2x2(A, B, X);
      real_t dydx = X[0], d2ydx2 = X[1];

      table->fit_data[0] = dydx;
      table->fit_data[1] = d2ydx2;
    }

    // Interior points are just quadratic fits about the second of 3 points.
    for (int i = 1; i < num_values-2; ++i)
    {
      real_t y1 = table->values[i-1];
      real_t y2 = table->values[i];
      real_t y3 = table->values[i+1];

      real_t A[4] = {-dx, dx, dx*dx, dx*dx};
      real_t B[2] = {y1-y2, y3-y2};
      real_t X[2];
      solve_2x2(A, B, X);
      real_t dydx = X[0], d2ydx2 = X[1];

      table->fit_data[2*i]   = dydx;
      table->fit_data[2*i+1] = d2ydx2;
    }

    // The last interval is a quadratic fit about the last of 3 points.
    {
      real_t y1 = table->values[num_values-3];
      real_t y2 = table->values[num_values-2];
      real_t y3 = table->values[num_values-1];

      real_t A[4] = {-dx, -2*dx, -dx*dx, -4*dx*dx};
      real_t B[2] = {y2-y3, y1-y3};
      real_t X[2];
      solve_2x2(A, B, X);
      real_t dydx = X[0], d2ydx2 = X[1];

      table->fit_data[2*(num_values-2)]   = dydx;
      table->fit_data[2*(num_values-2)+1] = d2ydx2;
    }

    // Padding.
    table->fit_data[2*num_values-2] = 0.0;
    table->fit_data[2*num_values-1] = 0.0;
  }
  return table;
}

void lookup1_free(lookup1_t* table)
{
  polymec_free(table->fit_data);
  polymec_free(table->values);
  polymec_free(table);
}

int lookup1_num_values(lookup1_t* table)
{
  return table->num_values;
}

lookup1_interp_t lookup1_interpolation(lookup1_t* table)
{
  return table->interp;
}

void lookup1_get_bounds(lookup1_t* table, real_t* x_min, real_t* x_max)
{
  *x_min = table->x_min;
  *x_max = table->x_max;
}

void lookup1_get_extreme_values(lookup1_t* table, real_t* y_min, real_t* y_max)
{
  *y_min = table->values[0];
  *y_max = table->values[table->num_values-1];
}

real_t lookup1_value(lookup1_t* table, real_t x)
{
  return table->compute(table, x);
}

