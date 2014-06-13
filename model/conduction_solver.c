// Copyright (c) 2012-2014, Jeffrey N. Johnson
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

#include "model/conduction_solver.h"
#include "core/unordered_map.h"

typedef struct
{
  mesh_t* mesh;
  real_t *lambda, *f;
  string_ptr_unordered_map_t* bcs;
  int_int_unordered_map_t* bcell_map;
} cond_t;

typedef struct
{
  real_t *alpha, *beta, *gamma;
  int N;
} cond_bc_t;

static cond_bc_t* cond_bc_new(real_t* alpha, real_t* beta, real_t* gamma, int N)
{
  ASSERT(N > 0);
  ASSERT((alpha != NULL) || (beta != NULL) || (gamma != NULL));

  cond_bc_t* bc = polymec_malloc(sizeof(cond_bc_t));
  bc->alpha = alpha;
  bc->beta = beta;
  bc->gamma = gamma;
  bc->N = N;
  return bc;
}

static void cond_bc_free(cond_bc_t* bc)
{
  polymec_free(bc);
}

static void kv_dtor(char* key, void* val)
{
  polymec_free(key);
  cond_bc_t* bc = val;
  cond_bc_free(bc);
}

// The dot product of u o T o v, where u and v are vectors and T is a 
// symmetric tensor.
static real_t dot(vector_t* u, real_t* T, vector_t* v)
{
  vector_t T_o_v = {.x = T[0] * v->x + T[1] * v->y + T[2] * v->z,
                    .y = T[1] * v->x + T[3] * v->y + T[4] * v->z,
                    .z = T[2] * v->x + T[4] * v->y + T[5] * v->z};
  return vector_dot(u, &T_o_v);
}

static int cond_ax(void* context, real_t* x, real_t* Ax, int N)
{
  // x and Ax are both cell centered. We now form the product Ax on the 
  // cell centers.

  // FIXME: Parallel exchange of fields happens here.

  cond_t* cond = context;
  mesh_t* mesh = cond->mesh;
  ASSERT(mesh->num_cells + cond->bcell_map->size == N);
  for (int cell = 0; cell < mesh->num_cells; ++cell)
  {
    int pos = 0, face;

    // Evaluate div ( lambda o grad phi) using the discrete divergence theorem.
    real_t div_lambda_grad = 0.0;
    point_t* xc = &mesh->cell_centers[cell];
    while (mesh_cell_next_face(mesh, cell, &pos, &face))
    {
      int neighbor = mesh_face_opp_cell(mesh, face, cell);
      real_t grad;
      real_t lambda = 1.0;
      if (neighbor != -1)
      {
        point_t* xn = &mesh->cell_centers[neighbor];
        grad = (x[neighbor] - x[cell]) / point_distance(xn, xc);
        if (cond->lambda != NULL)
          lambda = 0.5 * (cond->lambda[cell] + cond->lambda[neighbor]);
      }
      else
      {
        int bcell = *int_int_unordered_map_get(cond->bcell_map, face);
        point_t* xf = &mesh->face_centers[face];
        grad = (x[bcell] - x[cell]) / (2.0 * point_distance(xf, xc));
        if (cond->lambda != NULL)
          lambda = cond->lambda[cell];
      }
      real_t area = mesh->face_areas[face];
      div_lambda_grad += area * lambda * grad;
    }
    Ax[cell] = div_lambda_grad;
  }

  // Handle boundary conditions.
  int pos = 0;
  char* face_tag;
  void* bc_ptr;
  while (string_ptr_unordered_map_next(cond->bcs, &pos, &face_tag, &bc_ptr))
  {
    cond_bc_t* bc = bc_ptr;
    int num_faces;
    int* tag = mesh_tag(cond->mesh->face_tags, face_tag, &num_faces);
    for (int f = 0; f < num_faces; ++f)
    {
      int face = tag[f];
      int cell = mesh->face_cells[2*face];
      int bcell = *int_int_unordered_map_get(cond->bcell_map, face);

      // Construct this face's contribution to the LHS.
      real_t alpha = (bc->alpha != NULL) ? bc->alpha[f] : 0.0;
      real_t beta = (bc->beta != NULL) ? bc->beta[f] : 0.0;
      point_t* xc = &mesh->cell_centers[cell];
      point_t* xf = &mesh->face_centers[face];
      real_t grad = (x[bcell] - x[cell]) / (2.0 * point_distance(xc, xf));
      Ax[bcell] = alpha * 0.5 * (x[cell] + x[bcell]) + beta * grad;
    }
  }

  return 0;
}

static void cond_b(void* context, real_t* B, int N)
{
  cond_t* cond = context;
  mesh_t* mesh = cond->mesh;
  if (cond->f != NULL)
  {
    ASSERT(mesh->num_cells == N);
    for (int c = 0; c < N; ++c)
      B[c] = mesh->cell_volumes[c] * cond->f[c];
  }
  else
    memset(B, 0, sizeof(real_t) * N);

  // Handle boundary conditions.
  int pos = 0;
  char* face_tag;
  void* bc_ptr;
  while (string_ptr_unordered_map_next(cond->bcs, &pos, &face_tag, &bc_ptr))
  {
    cond_bc_t* bc = bc_ptr;
    int num_faces;
    int* tag = mesh_tag(cond->mesh->face_tags, face_tag, &num_faces);
    for (int f = 0; f < num_faces; ++f)
    {
      int face = tag[f];
      int bcell = *int_int_unordered_map_get(cond->bcell_map, face);

      // Construct this face's contribution to the RHS.
      real_t gamma = (bc->gamma != NULL) ? bc->gamma[f] : 0.0;
      B[bcell] = gamma;
    }
  }
}

static void cond_dtor(void* context)
{
  cond_t* cond = context;
  int_int_unordered_map_free(cond->bcell_map);
  string_ptr_unordered_map_free(cond->bcs);
  polymec_free(cond);
}

// This helper generates a mapping of faces to indices for imaginary "boundary" 
// cells used to enforce boundary conditions. These boundary cells are not to 
// be confused with ghost cells, which only exist for dealing with parallel 
// process boundaries.
static int_int_unordered_map_t* generate_boundary_cell_map(mesh_t* mesh)
{
  int_int_unordered_map_t* bcell_map = int_int_unordered_map_new();
  int index = mesh->num_cells;
  for (int cell = 0; cell < mesh->num_cells; ++cell)
  {
    int pos = 0, face;
    while (mesh_cell_next_face(mesh, cell, &pos, &face))
    {
      if (mesh_face_opp_cell(mesh, face, cell) == -1)
        int_int_unordered_map_insert(bcell_map, face, index++);
    }
  }
  return bcell_map;
}

krylov_solver_t* gmres_conduction_solver_new(mesh_t* mesh, real_t* lambda, real_t* f,
                                            int max_krylov_dim, int max_restarts)
{
  cond_t* cond = polymec_malloc(sizeof(cond_t));
  cond->mesh = mesh;
  cond->lambda = lambda;
  cond->f = f;
  cond->bcs = string_ptr_unordered_map_new();
  cond->bcell_map = generate_boundary_cell_map(mesh);
  int N = mesh->num_cells + cond->bcell_map->size;
  krylov_solver_vtable vtable = {.ax = cond_ax, .b = cond_b, .dtor = cond_dtor};
  return gmres_krylov_solver_new(mesh->comm, cond, vtable, N, max_krylov_dim, max_restarts);
}

krylov_solver_t* bicgstab_conduction_solver_new(mesh_t* mesh, real_t* lambda, real_t* f,
                                               int max_krylov_dim)
{
  cond_t* cond = polymec_malloc(sizeof(cond_t));
  cond->mesh = mesh;
  cond->lambda = lambda;
  cond->f = f;
  cond->bcs = string_ptr_unordered_map_new();
  cond->bcell_map = generate_boundary_cell_map(mesh);
  int N = mesh->num_cells + cond->bcell_map->size;
  krylov_solver_vtable vtable = {.ax = cond_ax, .b = cond_b, .dtor = cond_dtor};
  return bicgstab_krylov_solver_new(mesh->comm, cond, vtable, N, max_krylov_dim);
}

krylov_solver_t* tfqmr_conduction_solver_new(mesh_t* mesh, real_t* lambda, real_t* f,
                                            int max_krylov_dim)
{
  cond_t* cond = polymec_malloc(sizeof(cond_t));
  cond->mesh = mesh;
  cond->lambda = lambda;
  cond->f = f;
  cond->bcs = string_ptr_unordered_map_new();
  cond->bcell_map = generate_boundary_cell_map(mesh);
  int N = mesh->num_cells + cond->bcell_map->size;
  krylov_solver_vtable vtable = {.ax = cond_ax, .b = cond_b, .dtor = cond_dtor};
  return tfqmr_krylov_solver_new(mesh->comm, cond, vtable, N, max_krylov_dim);
}

void conduction_solver_add_bc(krylov_solver_t* conduction, const char* face_tag,
                              real_t* alpha, real_t* beta, real_t* gamma)
{
  cond_t* cond = krylov_solver_context(conduction);

  // Retrieve the mesh's face tag.
  ASSERT(mesh_has_tag(cond->mesh->face_tags, face_tag));
  int num_faces;
  mesh_tag(cond->mesh->face_tags, face_tag, &num_faces);

  // Set up a boundary condition and stick it into our map of BCs.
  cond_bc_t* bc = cond_bc_new(alpha, beta, gamma, num_faces);
  char* tag_name = string_dup(face_tag);
  string_ptr_unordered_map_insert_with_kv_dtor(cond->bcs, tag_name, bc, kv_dtor);
}

