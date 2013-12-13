// Copyright (c) 2012-2013, Jeffrey N. Johnson
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this 
// list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice, 
// this list of conditions and the following disclaimer in the documentation 
// and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#include "core/constant_st_func.h"
#include "core/linear_algebra.h"
#include "geometry/plane.h"

typedef struct
{
  // Plane parameters.
  vector_t n;
  point_t x;

  // Basis vectors (for projections).
  vector_t e1, e2, e3;
} plane_t;

static void plane_eval(void* ctx, point_t* x, double* result)
{
  plane_t* p = ctx;
  vector_t D;
  point_displacement(x, &p->x, &D);
  result[0] = vector_dot(&p->n, &D);
}

sp_func_t* plane_new(vector_t* n, point_t* x)
{
  // Set up a plane signed distance function.
  plane_t* p = malloc(sizeof(plane_t));
  sp_vtable vtable = {.eval = plane_eval, .dtor = free};
  sp_func_t* plane = sp_func_new("Plane (uninitialized)", p, vtable, SP_INHOMOGENEOUS, 1);

  // Set up all the innards.
  plane_reset(plane, n, x);
  return plane;
}

sp_func_t* plane_new_from_points(point_t* p1, point_t* p2, point_t* p3)
{
  ASSERT(!points_are_colinear(p1, p2, p3));

  vector_t e1, e2, n;
  point_displacement(p1, p2, &e1);
  point_displacement(p1, p3, &e2);
  vector_cross(&e1, &e2, &n);
  vector_normalize(&n);
  return plane_new(&n, p1);
}

sp_func_t* plane_new_best_fit(point_t* points, int num_points)
{
  ASSERT(num_points >= 3);

  if (num_points == 3)
  {
    // If exactly 3 points are given, we can exactly solve the equation of 
    // the plane.
    return plane_new_from_points(&points[0], &points[1], &points[2]);
  }
  else
  {
    // Otherwise, we resort to an orthogonal distance regression.

    // Note: the plane point x0 is the centroid of the given points.
    point_t x0 = {0.0, 0.0, 0.0};
    for (int i = 0; i < num_points; ++i)
    {
      x0.x += points[i].x;
      x0.y += points[i].y;
      x0.z += points[i].z;
    }
    x0.x /= num_points;
    x0.y /= num_points;
    x0.z /= num_points;

    double M[num_points*3]; // num_points x 3 matrix in column-major form.
    for (int i = 0; i < num_points; ++i)
    {
      M[3*i]   = points[i].x - x0.x;
      M[3*i+1] = points[i].y - x0.y;
      M[3*i+2] = points[i].z - x0.z;
    }

    // Compute the singular value decomposition of M to find the singular 
    // vectors. We will destroy M in the process and only compute the right 
    // singular vectors.
    double sigma[3]; // Singular values.
    double VT[3*3]; // Right singular vectors.
    int work_size = 10*num_points;
    double work[work_size]; // Work array.
    int one = 1, three = 3, info;
    char jobU = 'N', jobVT = 'S';
    dgesvd(&jobU, &jobVT, &num_points, &three, M, &num_points, sigma, 
           NULL, &one, VT, &three, work, &work_size, &info);
    ASSERT(info == 0);

    // The normal vector of the plane is the singular vector corresponding to 
    // the smallest singular value.
    double min_s = FLT_MAX;
    int min_index = -1;
    for (int i = 0; i < 3; ++i)
    {
      if (sigma[i] < min_s)
      {
        min_s = sigma[i];
        min_index = i;
      }
    }
    vector_t n;
    n.x = VT[3*min_index];
    n.y = VT[3*min_index+1];
    n.z = VT[3*min_index+2];
    return plane_new(&n, &x0);
  }
}

void plane_reset(sp_func_t* plane, vector_t* n, point_t* x)
{
  plane_t* p = sp_func_context(plane);
  vector_copy(&p->n, n);
  vector_normalize(&p->n);
  point_copy(&p->x, x);

  char plane_str[1024];
  sprintf(plane_str, "Plane (n = (%g, %g, %g), x = (%g, %g, %g))", 
          n->x, n->y, n->z, x->x, x->y, x->z);
  sp_func_rename(plane, plane_str);

  // Register the negative of the normal as the derivative.
  double nn[3];
  nn[0] = -n->x, nn[1] = -n->y, nn[2] = -n->z; 
  sp_func_t* G = constant_sp_func_new(3, nn);
  sp_func_register_deriv(plane, 1, G);

  // Set up our basis vectors.
  p->e3.x = nn[0], p->e3.y = nn[1], p->e3.z = nn[2];
  vector_normalize(&p->e3);
  compute_orthonormal_basis(&p->e3, &p->e1, &p->e2);
}

void plane_project(sp_func_t* plane, point_t* x, point2_t* xi)
{
  plane_t* p = sp_func_context(plane);
  vector_t v = {.x = x->x - p->x.x, .y = x->y - p->x.y, .z = x->z - p->x.z};
  double voe3 = vector_dot(&v, &p->e3);
  vector_t v_perp = {.x = v.x - voe3 * p->e3.x, .y = v.y - voe3 * p->e3.y, .z = v.z - voe3 * p->e3.z};
  xi->x = vector_dot(&v_perp, &p->e1);
  xi->y = vector_dot(&v_perp, &p->e2);
}

void plane_project2(sp_func_t* plane, point_t* x, double* eta, double* xi)
{
  plane_t* p = sp_func_context(plane);
  vector_t v = {.x = x->x - p->x.x, .y = x->y - p->x.y, .z = x->z - p->x.z};
  double voe3 = vector_dot(&v, &p->e3);
  vector_t v_perp = {.x = v.x - voe3 * p->e3.x, .y = v.y - voe3 * p->e3.y, .z = v.z - voe3 * p->e3.z};
  *eta = vector_dot(&v_perp, &p->e1);
  *xi = vector_dot(&v_perp, &p->e2);
}

void plane_embed(sp_func_t* plane, point2_t* xi, point_t* x)
{
  plane_t* p = sp_func_context(plane);
  x->x = p->x.x + p->e1.x * xi->x + p->e2.x * xi->y;
  x->y = p->x.y + p->e1.y * xi->x * p->e2.y + xi->y;
  x->z = p->x.z + p->e1.z * xi->x * p->e2.z + xi->y;
}

double plane_intersect_with_line(sp_func_t* plane, point_t* x0, vector_t* t)
{
  plane_t* p = sp_func_context(plane);
  double not = vector_dot(&p->n, t);
  if (not == 0.0) // No intersection!
    return -FLT_MAX;
  else
  {
    double numer = p->n.x * (p->x.x - x0->x) +
                   p->n.y * (p->x.y - x0->y) +
                   p->n.z * (p->x.z - x0->z);
    return numer / not;
  }
}

