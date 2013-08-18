// Copyright 2012-2013 Jeffrey Johnson.
// 
// This file is part of Polymec, and is licensed under the Apache License, 
// Version 2.0 (the "License"); you may not use this file except in 
// compliance with the License. You may may find the text of the license in 
// the LICENSE file at the top-level source directory, or obtain a copy of 
// it at
// 
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "integrators/rk_integrator.h"

typedef struct
{
  int order;
  void* context;
  rk_compute_deriv compute_deriv;
  integrator_dtor dtor;
  double *a, *b, *c; // Butcher tableau coefficients.
  double *y, **k;
} rk_t;

static void rk_dtor(void* context)
{
  rk_t* rk = context;
  if ((rk->context != NULL) && (rk->dtor != NULL))
    rk->dtor(rk->context);
  for (int s = 0; s < rk->order; ++s)
    free(rk->k[s]);
  free(rk->k);
  free(rk->y);
  free(rk);
}

static void rk_init(void* context, int N)
{
  rk_t* rk = context;
  rk->y = malloc(sizeof(double)*N);
  rk->k = malloc(sizeof(double*));
  for (int s = 0; s < rk->order; ++s) 
    rk->k[s] = malloc(sizeof(double)*N);
}

static void rk_step(void* context, double t1, double t2, double* solution, int N)
{
  ASSERT(t2 > t1);
  double h = t2 - t1;
  rk_t* rk = context;
  int order = rk->order;
  for (int i = 0; i < order; ++i)
  {
    double t = t1 + rk->c[i]*h;
    for (int n = 0; n < N; ++n)
      rk->y[n] = solution[n];
    for (int j = 0; j < i; ++j)
    {
      for (int n = 0; n < N; ++n)
        rk->y[n] += h * rk->a[order*i+j] * rk->k[j][n];
    }
    rk->compute_deriv(rk->context, t, rk->y, rk->k[i]);
  }

  for (int i = 0; i < order; ++i)
  {
    for (int n = 0; n < N; ++n)
      solution[n] += h * rk->b[i] * rk->k[i][n];
  }
}

// -----------------
// Butcher tableaus.
// -----------------

// 1st-order "Euler" method.
static double rk1_a[] = { 0.0 };
static double rk1_b[] = { 1.0 };
static double rk1_c[] = { 0.0 };

// 2nd-order Midpoint method.
static double rk2_a[] = { 0.0, 0.0, 
                          0.5, 0.0 };
static double rk2_b[] = { 0.0, 1.0 };
static double rk2_c[] = { 0.0, 0.5 };

// Kutta's 3rd-order method.
static double rk3_a[] = { 0.0, 0.0, 0.0, 
                          0.5, 0.0, 0.0, 
                         -1.0, 2.0, 0.0 };
static double rk3_b[] = { 1.0/6.0, 2.0/3.0, 1.0/6.0 };
static double rk3_c[] = { 0.0, 0.5, 1.0 };

// Classic 4th-order method.
static double rk4_a[] = { 0.0, 0.0, 0.0, 0.0, 
                          0.5, 0.0, 0.0, 0.0, 
                          0.0, 0.5, 0.0, 0.0,
                          0.0, 0.0, 1.0, 0.0 };
static double rk4_b[] = { 1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0 };
static double rk4_c[] = { 0.0, 0.5, 0.5, 1.0 };

integrator_t* rk_integrator_new(int order, 
                                void* context, 
                                rk_compute_deriv compute_deriv,
                                integrator_dtor dtor)
{
  ASSERT(order >= 1);
  ASSERT(order <= 4);
  ASSERT(compute_deriv != NULL);
  rk_t* rk = malloc(sizeof(rk_t));
  rk->order = order;
  rk->context = context;
  switch (order)
  {
    case 1: rk->a = rk1_a;
            rk->b = rk1_b;
            rk->c = rk1_c;
            break;
    case 2: rk->a = rk2_a;
            rk->b = rk2_b;
            rk->c = rk2_c;
            break;
    case 3: rk->a = rk3_a;
            rk->b = rk3_b;
            rk->c = rk3_c;
            break;
    case 4: rk->a = rk4_a;
            rk->b = rk4_b;
            rk->c = rk4_c;
            break;
    default: break;
  }
  rk->compute_deriv = compute_deriv;
  rk->dtor = dtor;
  integrator_vtable vtable = {.init = rk_init, .step = rk_step, .dtor = rk_dtor};
  char name[1024];
  sprintf(name, "RK%d", order);
  return integrator_new(name, rk, vtable, order, INTEGRATOR_EXPLICIT);
}

