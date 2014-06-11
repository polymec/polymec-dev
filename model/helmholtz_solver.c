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

static int helm_ax(void* context, real_t* x, real_t* Ax)
{
  // x and Ax are both cell centered. We now form the product Ax on the 
  // cell centers.

  // FIXME: Parallel exchange of fields happens here.

  helm_t* helm = context;
  mesh_t* mesh = helm->mesh;
  for (int cell = 0; cell < mesh->num_cells; ++cell)
  {
    int pos = 0, face;

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
        vector_t grad = {.x = dx / (xn->x - xc->x),            
                         .y = dx / (xn->y - xc->y),
                         .z = dx / (xn->z - xc->z)};
        div_grad += mesh->face_areas[face] * vector_dot(&mesh->face_normals[face], &grad);
      }
    }
    real_t k = helm->k[cell];
    Ax[cell] = div_grad + k * k * x[cell];
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

      // Construct this face's contribution to its cell's Laplacian.
      real_t alpha = (bc->alpha != NULL) ? bc->alpha[f] : 0.0;
      real_t beta = (bc->beta != NULL) ? bc->beta[f] : 0.0;
      real_t gamma = (bc->gamma != NULL) ? bc->gamma[f] : 0.0;
      real_t area = mesh->face_areas[face];
      vector_t* normal = &mesh->face_normals[face];
      int cell = mesh->face_cells[2*face];
      point_t* xc = &mesh->cell_centers[cell];
      point_t* xf = &mesh->face_centers[face];
      real_t phi_c = x[cell];
      real_t contrib;
      if (alpha == 0.0)
        contrib = area * gamma / beta;
      else 
      {
        // Compute phi_f, the solution on the face, from phi_c, the solution 
        // on the cell (and other things).
        real_t n_term = normal->x / (xf->x - xc->x) + 
                        normal->y / (xf->y - xc->y) + 
                        normal->z / (xc->z - xc->z);
        real_t phi_f = (beta * n_term * phi_c + gamma) / (alpha + beta * n_term);

        // Now construct the gradient and the contribution.
        vector_t grad = {.x = (phi_f - phi_c) / (xf->x - xc->x),            
                         .y = (phi_f - phi_c) / (xf->y - xc->y),
                         .z = (phi_f - phi_c) / (xf->z - xc->z)};
        contrib = area * vector_dot(normal, &grad);
      }
      Ax[cell] += contrib;
    }
  }

  return 0;
}

static void helm_free(void* context)
{
  helm_t* helm = context;
  string_ptr_unordered_map_free(helm->bcs);
  polymec_free(helm);
}

krylov_solver_t* helmholtz_solver_new(mesh_t* mesh, real_t* k, real_t* f)
{
  helm_t* helm = polymec_malloc(sizeof(helm_t));
  helm->mesh = mesh;
  helm->k = k;
  helm->f = f;
  helm->bcs = string_ptr_unordered_map_new();
  krylov_solver_vtable vtable = {.ax = helm_ax, .dtor = helm_free};

  // For now, we use GMRES.
  int max_krylov_dim = 15;
  int max_restarts = 5;
  return gmres_krylov_solver_new(mesh->comm, helm, vtable, mesh->num_cells,
                                 max_krylov_dim, max_restarts);
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

