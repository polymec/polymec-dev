// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "integrators/gauss_rules.h"
#include "meshless/mlpg_quadrature.h"

typedef struct
{
  point_cloud_t* cloud;
  real_t* extents;
  int N;
  real_t ratio; 
} mlpg_simple_t;

static mlpg_simple_t* mlpg_simple_new(point_cloud_t* cloud,
                                      real_t* extents,
                                      int num_points,
                                      real_t ratio)
{
  ASSERT(num_points > 0);
  ASSERT(ratio > 0.0);

  mlpg_simple_t* mlpg = polymec_malloc(sizeof(mlpg_simple_t));
  mlpg->cloud = cloud;
  mlpg->extents = extents;
  mlpg->N = num_points;
  mlpg->ratio = ratio;
  return mlpg;
}

static void mlpg_simple_free(void* context)
{
  polymec_free(context);
}

static int cube_surf_num_quad_points(void* context, int i)
{
  mlpg_simple_t* mlpg = context;
  return 6 * mlpg->N * mlpg->N;
}

static void cube_surf_get_quad(void* context, int i, point_t* points, real_t* weights, vector_t* normals)
{
  mlpg_simple_t* mlpg = context;
  real_t gauss_pts[mlpg->N], gauss_wts[mlpg->N];
  get_gauss_legendre_points(mlpg->N, gauss_wts, gauss_wts);
  point_t* xi = &mlpg->cloud->points[i];
  real_t hi = mlpg->extents[i];
  real_t L = mlpg->ratio * hi;
  real_t a = xi->x - 0.5*L, b = xi->x + 0.5*L;
  real_t c = xi->y - 0.5*L, d = xi->y + 0.5*L;
  real_t e = xi->z - 0.5*L, f = xi->z + 0.5*L;
  int N = mlpg->N;
  int m = 0;
  for (int jj = 0; jj < mlpg->N; ++jj)
  {
    for (int kk = 0; kk < mlpg->N; ++kk, ++m)
    {
      // X coordinates for each of the faces.
      points[m].x       = a;
      points[m+N*N].x   = b;
      real_t xj = (b - a) * gauss_pts[jj] + a;
      points[m+2*N*N].x = xj;
      points[m+3*N*N].x = xj;
      real_t xk = (b - a) * gauss_pts[kk] + a;
      points[m+4*N*N].x = xk;
      points[m+5*N*N].x = xk;

      // Y coordinates for each of the faces.
      real_t yj = (d - c) * gauss_pts[jj] + c;
      points[m].y = yj;
      points[m+N*N].y = yj;
      points[m+2*N*N].y = c;
      points[m+3*N*N].y = d;
      real_t yk = (f - e) * gauss_pts[kk] + c;
      points[m+4*N*N].y = yk;
      points[m+5*N*N].y = yk;

      // Z coordinates for each of the faces.
      real_t zj = (f - e) * gauss_pts[jj] + e;
      points[m].z = zj;
      points[m+N*N].z = zj;
      real_t zk = (f - e) * gauss_pts[kk] + e;
      points[m+2*N*N].z = zk;
      points[m+3*N*N].z = zk;
      points[m+4*N*N].y = e;
      points[m+5*N*N].y = f;
    }
  }
}

surface_integral_t* mlpg_cube_surface_integral_new(point_cloud_t* cloud,
                                                   real_t* extents,
                                                   int num_points,
                                                   real_t side_to_extent_ratio)
{
  mlpg_simple_t* mlpg = mlpg_simple_new(cloud, extents, num_points, 
                                        side_to_extent_ratio);
  surface_integral_vtable vtable = {.num_quad_points = cube_surf_num_quad_points,
                                    .get_quadrature = cube_surf_get_quad,
                                    .dtor = mlpg_simple_free};
  char name[1025];
  snprintf(name, 1024, "MLPG cube surface integral (N = %d, side/extent = %g)", 
           num_points, side_to_extent_ratio);
  return surface_integral_new(name, mlpg, vtable);
}

static int cube_vol_num_quad_points(void* context, int i)
{
  mlpg_simple_t* mlpg = context;
  return mlpg->N * mlpg->N * mlpg->N;
}

static void cube_vol_get_quad(void* context, int i, point_t* points, real_t* weights)
{
  mlpg_simple_t* mlpg = context;
  real_t gauss_pts[mlpg->N], gauss_wts[mlpg->N];
  get_gauss_legendre_points(mlpg->N, gauss_wts, gauss_wts);
  point_t* xi = &mlpg->cloud->points[i];
  real_t hi = mlpg->extents[i];
  real_t L = mlpg->ratio * hi;
  real_t a = xi->x - 0.5*L, b = xi->x + 0.5*L;
  real_t c = xi->y - 0.5*L, d = xi->y + 0.5*L;
  real_t e = xi->z - 0.5*L, f = xi->z + 0.5*L;
  real_t V = (b - a) * (d - c) * (f - e);
  int m = 0;
  for (int ii = 0; ii < mlpg->N; ++ii)
  {
    for (int jj = 0; jj < mlpg->N; ++jj)
    {
      for (int kk = 0; kk < mlpg->N; ++kk, ++m)
      {
        points[m].x = (b - a) * gauss_pts[ii] + a;
        points[m].y = (d - c) * gauss_pts[jj] + c;
        points[m].z = (f - e) * gauss_pts[kk] + e;
        weights[m] = gauss_wts[ii]*gauss_wts[jj]*gauss_wts[kk] * V;
      }
    }
  }
}

volume_integral_t* mlpg_cube_volume_integral_new(point_cloud_t* cloud,
                                                 real_t* extents,
                                                 int num_points,
                                                 real_t side_to_extent_ratio)
{
  mlpg_simple_t* mlpg = mlpg_simple_new(cloud, extents, num_points, 
                                        side_to_extent_ratio);
  volume_integral_vtable vtable = {.num_quad_points = cube_vol_num_quad_points,
                                   .get_quadrature = cube_vol_get_quad,
                                   .dtor = mlpg_simple_free};
  char name[1025];
  snprintf(name, 1024, "MLPG cube volume integral (N = %d, side/extent = %g)", 
           num_points, side_to_extent_ratio);
  return volume_integral_new(name, mlpg, vtable);
}

static int sphere_surf_num_quad_points(void* context, int i)
{
  mlpg_simple_t* mlpg = context;
  return mlpg->N * mlpg->N;
}

static void sphere_surf_get_quad(void* context, int i, point_t* points, real_t* weights, vector_t* normals)
{
  mlpg_simple_t* mlpg = context;
  real_t gauss_pts[mlpg->N], gauss_wts[mlpg->N];
  get_gauss_legendre_points(mlpg->N, gauss_wts, gauss_wts);
  point_t* xi = &mlpg->cloud->points[i];
  real_t hi = mlpg->extents[i];
  real_t a = mlpg->ratio * hi;
  int m = 0;
  for (int jj = 0; jj < mlpg->N; ++jj)
  {
    real_t eta_j = gauss_pts[jj];
    real_t sin_pietaj = sin(M_PI * eta_j);
    for (int kk = 0; kk < mlpg->N; ++kk, ++m)
    {
      real_t gamma_k = gauss_pts[kk];
      points[m].x = xi->x + a * cos(M_PI * eta_j);
      points[m].y = xi->y + a * sin_pietaj * cos(2.0 * M_PI * gamma_k);
      points[m].z = xi->z + a * sin_pietaj * sin(2.0 * M_PI * gamma_k);
      weights[m] = 2.0*M_PI*M_PI*a*a*a*gauss_pts[i]*gauss_pts[i] * 
                   sin_pietaj * gauss_wts[jj]*gauss_wts[kk];
    }
  }
}

surface_integral_t* mlpg_sphere_surface_integral_new(point_cloud_t* cloud,
                                                     real_t* extents,
                                                     int num_points,
                                                     real_t radius_to_extent_ratio)
{
  mlpg_simple_t* mlpg = mlpg_simple_new(cloud, extents, num_points, 
                                        radius_to_extent_ratio);
  surface_integral_vtable vtable = {.num_quad_points = sphere_surf_num_quad_points,
                                    .get_quadrature = sphere_surf_get_quad,
                                    .dtor = mlpg_simple_free};
  char name[1025];
  snprintf(name, 1024, "MLPG sphere surface integral (N = %d, radius/extent = %g)", 
           num_points, radius_to_extent_ratio);
  return surface_integral_new(name, mlpg, vtable);
}

static int sphere_vol_num_quad_points(void* context, int i)
{
  mlpg_simple_t* mlpg = context;
  return mlpg->N * mlpg->N * mlpg->N;
}

static void sphere_vol_get_quad(void* context, int i, point_t* points, real_t* weights)
{
  mlpg_simple_t* mlpg = context;
  real_t gauss_pts[mlpg->N], gauss_wts[mlpg->N];
  get_gauss_legendre_points(mlpg->N, gauss_wts, gauss_wts);
  point_t* xi = &mlpg->cloud->points[i];
  real_t hi = mlpg->extents[i];
  real_t a = mlpg->ratio * hi;
  int m = 0;
  for (int ii = 0; ii < mlpg->N; ++ii)
  {
    real_t xi_i = gauss_pts[ii];
    for (int jj = 0; jj < mlpg->N; ++jj)
    {
      real_t eta_j = gauss_pts[jj];
      real_t sin_pietaj = sin(M_PI * eta_j);
      for (int kk = 0; kk < mlpg->N; ++kk, ++m)
      {
        real_t gamma_k = gauss_pts[kk];
        points[m].x = xi->x + a * xi_i * cos(M_PI * eta_j);
        points[m].y = xi->y + a * xi_i * sin_pietaj * cos(2.0 * M_PI * gamma_k);
        points[m].z = xi->z + a * xi_i * sin_pietaj * sin(2.0 * M_PI * gamma_k);
        weights[m] = 2.0*M_PI*M_PI*a*a*a*gauss_pts[i]*gauss_pts[i] * 
                     sin_pietaj * gauss_wts[ii]*gauss_wts[jj]*gauss_wts[kk];
      }
    }
  }
}

volume_integral_t* mlpg_sphere_volume_integral_new(point_cloud_t* cloud,
                                                   real_t* extents,
                                                   int num_points,
                                                   real_t radius_to_extent_ratio)
{
  mlpg_simple_t* mlpg = mlpg_simple_new(cloud, extents, num_points, 
                                        radius_to_extent_ratio);
  volume_integral_vtable vtable = {.num_quad_points = sphere_vol_num_quad_points,
                                   .get_quadrature = sphere_vol_get_quad,
                                   .dtor = mlpg_simple_free};
  char name[1025];
  snprintf(name, 1024, "MLPG sphere volume integral (N = %d, radius/extent = %g)", 
           num_points, radius_to_extent_ratio);
  return volume_integral_new(name, mlpg, vtable);
}
