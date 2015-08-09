// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

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
  // FIXME
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
  // FIXME
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
  // FIXME
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
  // FIXME
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
