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

#include "model/helmholtz_solver.h"

typedef struct
{
  mesh_t* mesh;
  real_t *k, *f;
  string_ptr_unordered_map_t* bcs;
} helm_t;

typedef struct
{
  real_t *alpha, *beta, *gamma;
  int N;
} helm_bc_t;

static helm_bc_t* helm_bc_new(real_t* alpha, real_t* beta, real_t* gamma, int N)
{
  ASSERT(N > 0);
  ASSERT((alpha != NULL) || (beta != NULL) || (gamma != NULL));

  helm_bc_t* bc = polymec_malloc(sizeof(helm_bc_t));
  bc->alpha = alpha;
  bc->beta = beta;
  bc->gamma = gamma;
  bc->N = N;
  return bc;
}

static void helm_bc_free(helm_bc_t* bc)
{
  polymec_free(bc);
}

static void kv_dtor(char* key, void* val)
{
  polymec_free(key);
  helm_bc_t* bc = val;
  helm_bc_free(bc);
}

static int helm_ax(void* context, real_t* x, real_t* Ax, int N)
{
  // x and Ax are both cell centered. We now form the product Ax on the 
  // cell centers.

  // FIXME: Parallel exchange of fields happens here.

  helm_t* helm = context;
  mesh_t* mesh = helm->mesh;
  ASSERT(mesh->num_cells == N);
  for (int cell = 0; cell < mesh->num_cells; ++cell)
  {
    int pos = 0, face;

    real_t volume = mesh->cell_volumes[cell];

    // Evaluate the Laplacian of phi using the discrete divergence theorem.
    real_t div_grad = 0.0;
    point_t* xc = &mesh->cell_centers[cell];
    while (mesh_cell_next_face(mesh, cell, &pos, &face))
    {
      int neighbor = mesh_face_opp_cell(mesh, face, cell);
      if (neighbor != -1)
      {
        real_t dx = x[neighbor] - x[cell];
        point_t* xn = &mesh->cell_centers[neighbor];
        real_t grad = dx / point_distance(xn, xc);           
        div_grad += mesh->face_areas[face] * grad;
      }
    }
    real_t k = (helm->k != NULL) ? helm->k[cell]: 0.0;
    Ax[cell] = div_grad + k * k * x[cell] * volume;
  }

  // Handle boundary conditions.
  int pos = 0;
  char* face_tag;
  void* bc_ptr;
  while (string_ptr_unordered_map_next(helm->bcs, &pos, &face_tag, &bc_ptr))
  {
    helm_bc_t* bc = bc_ptr;
    int num_faces;
    int* tag = mesh_tag(helm->mesh->face_tags, face_tag, &num_faces);
    for (int f = 0; f < num_faces; ++f)
    {
      int face = tag[f];

      // Construct this face's contribution to the LHS.
      real_t alpha = (bc->alpha != NULL) ? bc->alpha[f] : 0.0;
      real_t beta = (bc->beta != NULL) ? bc->beta[f] : 0.0;
      real_t gamma = (bc->gamma != NULL) ? bc->gamma[f] : 0.0;
      vector_t* normal = &mesh->face_normals[face];
      int cell = mesh->face_cells[2*face];
      point_t* xc = &mesh->cell_centers[cell];
      point_t* xf = &mesh->face_centers[face];
      real_t contrib;
      if (beta == 0.0)
      {
        // Supposing for the moment the existence of a ghost cell equidistant
        // from this face, we have 
        // alpha * phi_f = gamma,
        // phi_f = (phi_i + phi_g)/2.0,
        // phi_g = phi_i + 2 * (gamma - phi_i) / (xf - xc),
        // so alpha * (1 - 1/(xf-xc)) * phi_i = gamma * (1 - alpha / (xf - xc))
        contrib = alpha * (1.0 - 1.0 / point_distance(xf, xc));
      }
      else 
      {
#if 0
        // alpha * phi + beta * n o grad phi = gamma, and 
        // grad phi = (phi_f - phi_c) / (xf - xc)...
        real_t n_o_grad = normal->x / (xf->x - xc->x) + 
                          normal->y / (xf->y - xc->y) + 
                          normal->z / (xc->z - xc->z);
        // ...so phi_f = (beta * n_o_grad * phi_c + gamma) / (alpha + beta * n_o_grad).
        real_t phi_f_coeff = beta * n_o_grad / (alpha + beta * n_term);
        contrib = (1.0 + beta * n_o_grad * phi_f_coeff) * x[cell];
#endif
      }
      Ax[cell] += contrib;
    }
  }

  return 0;
}

static void helm_b(void* context, real_t* B, int N)
{
  helm_t* helm = context;
  mesh_t* mesh = helm->mesh;
  if (helm->f != NULL)
  {
    ASSERT(mesh->num_cells == N);
    for (int c = 0; c < N; ++c)
      B[c] = -mesh->cell_volumes[c] * helm->f[c];
  }
  else
    memset(B, 0, sizeof(real_t) * N);

  // Handle boundary conditions.
  int pos = 0;
  char* face_tag;
  void* bc_ptr;
  while (string_ptr_unordered_map_next(helm->bcs, &pos, &face_tag, &bc_ptr))
  {
    helm_bc_t* bc = bc_ptr;
    int num_faces;
    int* tag = mesh_tag(helm->mesh->face_tags, face_tag, &num_faces);
    for (int f = 0; f < num_faces; ++f)
    {
      int face = tag[f];

      // Construct this face's contribution to the RHS.
      real_t alpha = (bc->alpha != NULL) ? bc->alpha[f] : 0.0;
      real_t beta = (bc->beta != NULL) ? bc->beta[f] : 0.0;
      real_t gamma = (bc->gamma != NULL) ? bc->gamma[f] : 0.0;
      vector_t* normal = &mesh->face_normals[face];
      int cell = mesh->face_cells[2*face];
      point_t* xc = &mesh->cell_centers[cell];
      point_t* xf = &mesh->face_centers[face];
      real_t contrib;
      if (beta == 0.0)
      {
        // Supposing for the moment the existence of a ghost cell equidistant
        // from this face, we have 
        // alpha * phi_f = gamma,
        // phi_f = (phi_i + phi_g)/2.0,
        // phi_g = phi_i + 2 * (gamma - phi_i) / (xf - xc),
        // so alpha * (1 - 1/(xf-xc)) * phi_i = gamma * (1 - alpha / (xf - xc))
        contrib = gamma * (1.0 - alpha / point_distance(xf, xc));
      }
      else 
      {
#if 0
        // alpha * phi + beta * n o grad phi = gamma, and 
        // grad phi = (phi_f - phi_c) / (xf - xc)...
        real_t n_o_grad = normal->x / (xf->x - xc->x) + 
                          normal->y / (xf->y - xc->y) + 
                          normal->z / (xc->z - xc->z);
        // ...so phi_f = (beta * n_o_grad * phi_c + gamma) / (alpha + beta * n_o_grad).
        real_t phi_f_rhs = gamma / (alpha + beta * n_term);
        contrib = -beta * n_o_grad * phi_f_rhs;
#endif
      }
      B[cell] += contrib;
      printf("%d: alpha = %g, beta = %g, gamma = %g, B += %g\n", cell, alpha, beta, gamma, contrib);
    }
  }
}

static void helm_free(void* context)
{
  helm_t* helm = context;
  string_ptr_unordered_map_free(helm->bcs);
  polymec_free(helm);
}

krylov_solver_t* gmres_helmholtz_solver_new(mesh_t* mesh, real_t* k, real_t* f,
                                            int max_krylov_dim, int max_restarts)
{
  helm_t* helm = polymec_malloc(sizeof(helm_t));
  helm->mesh = mesh;
  helm->k = k;
  helm->f = f;
  helm->bcs = string_ptr_unordered_map_new();
  krylov_solver_vtable vtable = {.ax = helm_ax, .b = helm_b, .dtor = helm_free};
  return gmres_krylov_solver_new(mesh->comm, helm, vtable, mesh->num_cells,
                                 max_krylov_dim, max_restarts);
}

krylov_solver_t* bicgstab_helmholtz_solver_new(mesh_t* mesh, real_t* k, real_t* f,
                                               int max_krylov_dim)
{
  helm_t* helm = polymec_malloc(sizeof(helm_t));
  helm->mesh = mesh;
  helm->k = k;
  helm->f = f;
  helm->bcs = string_ptr_unordered_map_new();
  krylov_solver_vtable vtable = {.ax = helm_ax, .b = helm_b, .dtor = helm_free};
  return bicgstab_krylov_solver_new(mesh->comm, helm, vtable, mesh->num_cells,
                                    max_krylov_dim);
}

krylov_solver_t* tfqmr_helmholtz_solver_new(mesh_t* mesh, real_t* k, real_t* f,
                                            int max_krylov_dim)
{
  helm_t* helm = polymec_malloc(sizeof(helm_t));
  helm->mesh = mesh;
  helm->k = k;
  helm->f = f;
  helm->bcs = string_ptr_unordered_map_new();
  krylov_solver_vtable vtable = {.ax = helm_ax, .b = helm_b, .dtor = helm_free};
  return tfqmr_krylov_solver_new(mesh->comm, helm, vtable, mesh->num_cells,
                                 max_krylov_dim);
}

void helmholtz_solver_add_bc(krylov_solver_t* helmholtz, const char* face_tag,
                             real_t* alpha, real_t* beta, real_t* gamma)
{
  helm_t* helm = krylov_solver_context(helmholtz);

  // Retrieve the mesh's face tag.
  ASSERT(mesh_has_tag(helm->mesh->face_tags, face_tag));
  int num_faces;
  mesh_tag(helm->mesh->face_tags, face_tag, &num_faces);

  // Set up a boundary condition and stick it into our map of BCs.
  helm_bc_t* bc = helm_bc_new(alpha, beta, gamma, num_faces);
  char* tag_name = string_dup(face_tag);
  string_ptr_unordered_map_insert_with_kv_dtor(helm->bcs, tag_name, bc, kv_dtor);
}

