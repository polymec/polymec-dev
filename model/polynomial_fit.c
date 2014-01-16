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

#include "core/polynomial.h"
#include "core/linear_algebra.h"
#include "model/polynomial_fit.h"

struct polynomial_fit_t 
{
  char* name;
  void* context;
  polynomial_fit_vtable vtable;
  int num_comps;

  int degree;
  int point_index;
  int num_points;

  // Polynomials for fit.
  polynomial_t** polys;
};

polynomial_fit_t* polynomial_fit_new(const char* name,
                                     void* context,
                                     polynomial_fit_vtable vtable,
                                     int num_comps)
{
  ASSERT(num_comps > 0);
  ASSERT(vtable.num_interior_neighbors != NULL);
  ASSERT(vtable.get_interior_neighbors != NULL);
  ASSERT(vtable.get_interior_data != NULL);
  ASSERT(vtable.num_boundary_neighbors != NULL);
  ASSERT(vtable.get_boundary_neighbors != NULL);
  ASSERT(vtable.get_boundary_data != NULL);
  ASSERT(vtable.targeted_degree != NULL);

  polynomial_fit_t* fit = malloc(sizeof(polynomial_fit_t));
  fit->name = string_dup(name);
  fit->context = context;
  fit->vtable = vtable;
  fit->num_comps = num_comps;

  fit->degree = -1;
  fit->point_index = -1;
  fit->num_points = 0;
  fit->polys = NULL;

  return fit;
}

void polynomial_fit_free(polynomial_fit_t* fit)
{
  if (fit->polys != NULL)
  {
    for (int i = 0; i < fit->num_comps; ++i)
      fit->polys[i] = NULL;
  }
  if ((fit->context != NULL) && (fit->vtable.dtor != NULL))
    fit->vtable.dtor(fit->context);
  if (fit->name != NULL)
    free(fit->name);
  free(fit);
}

int polynomial_fit_num_comps(polynomial_fit_t* fit)
{
  return fit->num_comps;
}

void polynomial_fit_compute(polynomial_fit_t* fit, int point_index)
{
  ASSERT(point_index >= 0);
  fit->point_index = point_index;

  // Find the number of neighbors near the given point.
  int num_int_points = 1 + fit->vtable.num_interior_neighbors(fit->context, fit->point_index);
  int num_bnd_points = fit->vtable.num_boundary_neighbors(fit->context, fit->point_index);
  fit->num_points = num_int_points + num_bnd_points;

  // Determine the targeted degree for the number of point
  fit->degree = fit->vtable.targeted_degree(fit->context, fit->num_points);
  int dim = polynomial_basis_dim(fit->degree);
  ASSERT(dim >= fit->num_neighbors);

  // Get the indices of all the points.
  int int_indices[num_int_points];
  int_indices[0] = fit->point_index;
  fit->vtable.get_interior_neighbors(fit->context, fit->point_index, &int_indices[1]);
  int bnd_indices[num_bnd_points];
  fit->vtable.get_boundary_neighbors(fit->context, fit->point_index, bnd_indices);

  // Now fit each component.
  point_t int_points[num_int_points], bnd_points[num_bnd_points];
  real_t int_values[num_int_points], bnd_values[num_bnd_points];
  vector_t bnd_normals[num_bnd_points];
  real_t poly_coeffs[dim];
  for (int c = 0; c < fit->num_comps; ++c)
  {

    // Fetch the interior points and values.
    fit->vtable.get_interior_data(fit->context, c, int_indices, num_int_points,
                                  int_points, int_values);

    // Fetch the boundary points, normal vectors, and values.
    fit->vtable.get_boundary_data(fit->context, c, bnd_indices, num_bnd_points,
                                  bnd_points, bnd_normals, bnd_values);

    // Now fit the component data to a polynomial.
    fit->vtable.fit_component(fit->context, fit->degree, 
                              int_points, int_values, num_int_points,
                              bnd_points, bnd_normals, bnd_values, num_bnd_points,
                              poly_coeffs);

    // Construct the polynomial for this component.
    if (fit->polys[c] == NULL)
      fit->polys[c] = polynomial_new(fit->degree, poly_coeffs, &int_points[0]);
    else
    {
      memcpy(polynomial_coeffs(fit->polys[c]), poly_coeffs, sizeof(real_t) * dim);
      *polynomial_x0(fit->polys[c]) = int_points[0];
    }
  }
}

int polynomial_fit_degree(polynomial_fit_t* fit)
{
  return fit->degree;
}

void polynomial_fit_eval(polynomial_fit_t* fit, point_t* x, real_t* value)
{
  for (int i = 0; i < fit->num_comps; ++i)
    value[i] = polynomial_value(fit->polys[i], x);
}

void polynomial_fit_eval_deriv(polynomial_fit_t* fit, 
                               point_t* x, 
                               int x_deriv,
                               int y_deriv,
                               int z_deriv,
                               real_t* deriv)
{
  for (int i = 0; i < fit->num_comps; ++i)
    deriv[i] = polynomial_deriv_value(fit->polys[i], x_deriv, y_deriv, z_deriv, x);
}

//------------------------------------------------------------------------
//                Freebie polynomial fit constructors
//------------------------------------------------------------------------

// Cell-centered, fixed polynomial degree.
typedef struct 
{
  mesh_t* mesh;
  real_t* data;
  int degree;

  // Machinery for fit_component().
  void (*fit_component)(void*, int, point_t*, real_t*, int, point_t*, vector_t*, real_t*, int, real_t*);
  void* fit_component_context;
  void (*dtor)(void*);
} cc_fixed_degree_t;

static int cc_fixed_degree_num_interior_neighbors(void* context, int point_index)
{
  // FIXME
  return 0;
}

static void cc_fixed_degree_get_interior_neighbors(void* context, int point_index, int* neighbor_indices)
{
  // FIXME
}

static int cc_fixed_degree_num_boundary_neighbors(void* context, int point_index)
{
  // FIXME
  return 0;
}

static void cc_fixed_degree_get_boundary_neighbors(void* context, int point_index, int* neighbor_indices)
{
  // FIXME
}

static void cc_fixed_degree_get_interior_data(void* context, int component, int* point_indices, int num_points,
                                              point_t* point_coords, real_t* point_values)
{
  // FIXME
}

static void cc_fixed_degree_get_boundary_data(void* context, int component, int* point_indices, int num_points,
                                              point_t* point_coords, vector_t* boundary_normals, real_t* point_values)
{
  // FIXME
}

static int cc_fixed_degree_targeted_degree(void* context, int point_index)
{
  // FIXME
  return 0;
}

static void cc_fixed_degree_fit_component(void* context, int degree,
                                          point_t* interior_points, real_t* interior_values, int num_interior_points,
                                          point_t* boundary_points, vector_t* boundary_normals, real_t* boundary_values, int num_boundary_points,
                                          real_t* poly_coeffs)
{
  cc_fixed_degree_t* fit = context;
  fit->fit_component(fit->fit_component_context, degree, interior_points, interior_values, num_interior_points,
                     boundary_points, boundary_normals, boundary_values, num_boundary_points, poly_coeffs);
}

static void cc_fixed_degree_dtor(void* context)
{
  cc_fixed_degree_t* fit = context;
  if ((fit->fit_component_context != NULL) && (fit->dtor != NULL))
    fit->dtor(fit->fit_component_context);
}

polynomial_fit_t* cc_fixed_degree_polynomial_fit_new(int num_comps,
                                                     mesh_t* mesh,
                                                     real_t* data,
                                                     int degree,
                                                     void (*fit_component)(void*, int, point_t*, real_t*, int, point_t*, vector_t*, real_t*, int, real_t*),
                                                     void* fit_component_context,
                                                     void (*dtor)(void*))
{
  ASSERT(mesh != NULL);
  ASSERT(data != NULL);
  ASSERT(degree >= 0);

  char name[1024];
  snprintf(name, 1024, "cell-centered fixed degree (%d)", degree);
  cc_fixed_degree_t* context = malloc(sizeof(cc_fixed_degree_t));
  context->mesh = mesh;
  context->data = data;
  context->degree = degree;
  polynomial_fit_vtable vtable = {.num_interior_neighbors = cc_fixed_degree_num_interior_neighbors,
                                  .get_interior_neighbors = cc_fixed_degree_get_interior_neighbors,
                                  .get_interior_data = cc_fixed_degree_get_interior_data,
                                  .num_boundary_neighbors = cc_fixed_degree_num_boundary_neighbors,
                                  .get_boundary_neighbors = cc_fixed_degree_get_boundary_neighbors,
                                  .get_boundary_data = cc_fixed_degree_get_boundary_data,
                                  .targeted_degree = cc_fixed_degree_targeted_degree,
                                  .fit_component = cc_fixed_degree_fit_component,
                                  .dtor = cc_fixed_degree_dtor};
  return polynomial_fit_new(name, context, vtable, num_comps);
}

// Cell-centered, variable polynomial degree.
typedef struct 
{
  mesh_t* mesh;
  real_t* data;
  int depth;

  // Machinery for fit_component().
  void (*fit_component)(void*, int, point_t*, real_t*, int, point_t*, vector_t*, real_t*, int, real_t*);
  void* fit_component_context;
  void (*dtor)(void*);
} cc_var_degree_t;

static int cc_var_degree_num_interior_neighbors(void* context, int point_index)
{
  // FIXME
  return 0;
}

static void cc_var_degree_get_interior_neighbors(void* context, int point_index, int* neighbor_indices)
{
  // FIXME
}

static int cc_var_degree_num_boundary_neighbors(void* context, int point_index)
{
  // FIXME
  return 0;
}

static void cc_var_degree_get_boundary_neighbors(void* context, int point_index, int* neighbor_indices)
{
  // FIXME
}

static void cc_var_degree_get_interior_data(void* context, int component, int* point_indices, int num_points,
                                              point_t* point_coords, real_t* point_values)
{
  // FIXME
}

static void cc_var_degree_get_boundary_data(void* context, int component, int* point_indices, int num_points,
                                              point_t* point_coords, vector_t* boundary_normals, real_t* point_values)
{
  // FIXME
}

static int cc_var_degree_targeted_degree(void* context, int point_index)
{
  // FIXME
  return 0;
}

static void cc_var_degree_fit_component(void* context, int degree,
                                        point_t* interior_points, real_t* interior_values, int num_interior_points,
                                        point_t* boundary_points, vector_t* boundary_normals, real_t* boundary_values, int num_boundary_points,
                                        real_t* poly_coeffs)
{
  cc_var_degree_t* fit = context;
  fit->fit_component(fit->fit_component_context, degree, interior_points, interior_values, num_interior_points,
                     boundary_points, boundary_normals, boundary_values, num_boundary_points, poly_coeffs);
}

static void cc_var_degree_dtor(void* context)
{
  cc_var_degree_t* fit = context;
  if ((fit->fit_component_context != NULL) && (fit->dtor != NULL))
    fit->dtor(fit->fit_component_context);
}

polynomial_fit_t* cc_variable_degree_polynomial_fit_new(int num_comps,
                                                        mesh_t* mesh,
                                                        real_t* data,
                                                        int neighbor_search_depth,
                                                        void (*fit_component)(void*, int, point_t*, real_t*, int, point_t*, vector_t*, real_t*, int, real_t*),
                                                        void* fit_component_context,
                                                        void (*dtor)(void*))
{
  ASSERT(mesh != NULL);
  ASSERT(data != NULL);
  ASSERT(neighbor_search_depth >= 1);

  char name[1024];
  snprintf(name, 1024, "cell-centered variable degree (depth %d)", neighbor_search_depth);
  cc_var_degree_t* context = malloc(sizeof(cc_var_degree_t));
  context->mesh = mesh;
  context->data = data;
  context->depth = neighbor_search_depth;
  polynomial_fit_vtable vtable = {.num_interior_neighbors = cc_var_degree_num_interior_neighbors,
                                  .get_interior_neighbors = cc_var_degree_get_interior_neighbors,
                                  .get_interior_data = cc_var_degree_get_interior_data,
                                  .num_boundary_neighbors = cc_var_degree_num_boundary_neighbors,
                                  .get_boundary_neighbors = cc_var_degree_get_boundary_neighbors,
                                  .get_boundary_data = cc_var_degree_get_boundary_data,
                                  .targeted_degree = cc_var_degree_targeted_degree,
                                  .fit_component = cc_var_degree_fit_component,
                                  .dtor = cc_var_degree_dtor};
  return polynomial_fit_new(name, context, vtable, num_comps);
}

// Point cloud, fixed polynomial degree.
typedef struct 
{
  point_cloud_t* points;
  real_t* data;
  int degree;

  // Machinery for fit_component().
  void (*fit_component)(void*, int, point_t*, real_t*, int, point_t*, vector_t*, real_t*, int, real_t*);
  void* fit_component_context;
  void (*dtor)(void*);
} cloud_fixed_degree_t;

static int cloud_fixed_degree_num_interior_neighbors(void* context, int point_index)
{
  // FIXME
  return 0;
}

static void cloud_fixed_degree_get_interior_neighbors(void* context, int point_index, int* neighbor_indices)
{
  // FIXME
}

static int cloud_fixed_degree_num_boundary_neighbors(void* context, int point_index)
{
  // FIXME
  return 0;
}

static void cloud_fixed_degree_get_boundary_neighbors(void* context, int point_index, int* neighbor_indices)
{
  // FIXME
}

static void cloud_fixed_degree_get_interior_data(void* context, int component, int* point_indices, int num_points,
                                              point_t* point_coords, real_t* point_values)
{
  // FIXME
}

static void cloud_fixed_degree_get_boundary_data(void* context, int component, int* point_indices, int num_points,
                                              point_t* point_coords, vector_t* boundary_normals, real_t* point_values)
{
  // FIXME
}

static int cloud_fixed_degree_targeted_degree(void* context, int point_index)
{
  // FIXME
  return 0;
}
static void cloud_fixed_degree_fit_component(void* context, int degree,
                                             point_t* interior_points, real_t* interior_values, int num_interior_points,
                                             point_t* boundary_points, vector_t* boundary_normals, real_t* boundary_values, int num_boundary_points,
                                             real_t* poly_coeffs)
{
  cloud_fixed_degree_t* fit = context;
  fit->fit_component(fit->fit_component_context, degree, interior_points, interior_values, num_interior_points,
                     boundary_points, boundary_normals, boundary_values, num_boundary_points, poly_coeffs);
}

static void cloud_fixed_degree_dtor(void* context)
{
  cloud_fixed_degree_t* fit = context;
  if ((fit->fit_component_context != NULL) && (fit->dtor != NULL))
    fit->dtor(fit->fit_component_context);
}

polynomial_fit_t* point_cloud_fixed_degree_polynomial_fit_new(int num_comps,
                                                              point_cloud_t* points,
                                                              point_cloud_neighbor_search_t* search,
                                                              real_t* data,
                                                              int degree,
                                                              void (*fit_component)(void*, int, point_t*, real_t*, int, point_t*, vector_t*, real_t*, int, real_t*),
                                                              void* fit_component_context,
                                                              void (*dtor)(void*))
{
  ASSERT(points != NULL);
  ASSERT(search != NULL);
  ASSERT(data != NULL);
  ASSERT(degree >= 0);

  char name[1024];
  snprintf(name, 1024, "point cloud fixed degree (%d)", degree);
  cloud_fixed_degree_t* context = malloc(sizeof(cloud_fixed_degree_t));
  context->points = points;
  context->data = data;
  context->degree = degree;
  polynomial_fit_vtable vtable = {.num_interior_neighbors = cloud_fixed_degree_num_interior_neighbors,
                                  .get_interior_neighbors = cloud_fixed_degree_get_interior_neighbors,
                                  .get_interior_data = cloud_fixed_degree_get_interior_data,
                                  .num_boundary_neighbors = cloud_fixed_degree_num_boundary_neighbors,
                                  .get_boundary_neighbors = cloud_fixed_degree_get_boundary_neighbors,
                                  .get_boundary_data = cloud_fixed_degree_get_boundary_data,
                                  .targeted_degree = cloud_fixed_degree_targeted_degree,
                                  .fit_component = cloud_fixed_degree_fit_component,
                                  .dtor = cloud_fixed_degree_dtor};
  return polynomial_fit_new(name, context, vtable, num_comps);
}

// Point cloud, variable polynomial degree.
typedef struct 
{
  point_cloud_t* points;
  real_t* data;
  int depth;

  // Machinery for fit_component().
  void (*fit_component)(void*, int, point_t*, real_t*, int, point_t*, vector_t*, real_t*, int, real_t*);
  void* fit_component_context;
  void (*dtor)(void*);
} cloud_var_degree_t;

static int cloud_var_degree_num_interior_neighbors(void* context, int point_index)
{
  // FIXME
  return 0;
}

static void cloud_var_degree_get_interior_neighbors(void* context, int point_index, int* neighbor_indices)
{
  // FIXME
}

static int cloud_var_degree_num_boundary_neighbors(void* context, int point_index)
{
  // FIXME
  return 0;
}

static void cloud_var_degree_get_boundary_neighbors(void* context, int point_index, int* neighbor_indices)
{
  // FIXME
}

static void cloud_var_degree_get_interior_data(void* context, int component, int* point_indices, int num_points,
                                              point_t* point_coords, real_t* point_values)
{
  // FIXME
}

static void cloud_var_degree_get_boundary_data(void* context, int component, int* point_indices, int num_points,
                                              point_t* point_coords, vector_t* boundary_normals, real_t* point_values)
{
  // FIXME
}

static int cloud_var_degree_targeted_degree(void* context, int point_index)
{
  // FIXME
  return 0;
}

static void cloud_var_degree_fit_component(void* context, int degree,
                                           point_t* interior_points, real_t* interior_values, int num_interior_points,
                                           point_t* boundary_points, vector_t* boundary_normals, real_t* boundary_values, int num_boundary_points,
                                           real_t* poly_coeffs)
{
  cloud_var_degree_t* fit = context;
  fit->fit_component(fit->fit_component_context, degree, interior_points, interior_values, num_interior_points,
                     boundary_points, boundary_normals, boundary_values, num_boundary_points, poly_coeffs);
}

static void cloud_var_degree_dtor(void* context)
{
  cloud_var_degree_t* fit = context;
  if ((fit->fit_component_context != NULL) && (fit->dtor != NULL))
    fit->dtor(fit->fit_component_context);
}

polynomial_fit_t* point_cloud_variable_degree_polynomial_fit_new(int num_comps,
                                                                 point_cloud_t* points,
                                                                 point_cloud_neighbor_search_t* search,
                                                                 real_t* data,
                                                                 int neighbor_search_depth,
                                                                 void (*fit_component)(void*, int, point_t*, real_t*, int, point_t*, vector_t*, real_t*, int, real_t*),
                                                                 void* fit_component_context,
                                                                 void (*dtor)(void*))
{
  ASSERT(points != NULL);
  ASSERT(search != NULL);
  ASSERT(data != NULL);
  ASSERT(neighbor_search_depth >= 1);

  char name[1024];
  snprintf(name, 1024, "point cloud variable degree (depth = %d)", neighbor_search_depth);
  cloud_var_degree_t* context = malloc(sizeof(cloud_var_degree_t));
  context->points = points;
  context->data = data;
  context->depth = neighbor_search_depth;
  polynomial_fit_vtable vtable = {.num_interior_neighbors = cloud_var_degree_num_interior_neighbors,
                                  .get_interior_neighbors = cloud_var_degree_get_interior_neighbors,
                                  .get_interior_data = cloud_var_degree_get_interior_data,
                                  .num_boundary_neighbors = cloud_var_degree_num_boundary_neighbors,
                                  .get_boundary_neighbors = cloud_var_degree_get_boundary_neighbors,
                                  .get_boundary_data = cloud_var_degree_get_boundary_data,
                                  .targeted_degree = cloud_var_degree_targeted_degree,
                                  .fit_component = cloud_var_degree_fit_component,
                                  .dtor = cloud_var_degree_dtor};
  return polynomial_fit_new(name, context, vtable, num_comps);
}

