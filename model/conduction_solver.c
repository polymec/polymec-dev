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
#include "core/block_jacobi_preconditioner.h"

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
  int num_faces;
} cond_bc_t;

static cond_bc_t* cond_bc_new(int num_faces)
{
  ASSERT(num_faces > 0);

  cond_bc_t* bc = polymec_malloc(sizeof(cond_bc_t));
  bc->alpha = polymec_malloc(sizeof(real_t) * num_faces);
  memset(bc->alpha, 0, sizeof(real_t) * num_faces);
  bc->beta = polymec_malloc(sizeof(real_t) * num_faces);
  memset(bc->beta, 0, sizeof(real_t) * num_faces);
  bc->gamma = polymec_malloc(sizeof(real_t) * num_faces);
  memset(bc->gamma, 0, sizeof(real_t) * num_faces);
  bc->num_faces = num_faces;
  return bc;
}

static void cond_bc_free(cond_bc_t* bc)
{
  polymec_free(bc->alpha);
  polymec_free(bc->beta);
  polymec_free(bc->gamma);
  polymec_free(bc);
}

static void kv_dtor(char* key, void* val)
{
  polymec_free(key);
  cond_bc_t* bc = val;
  cond_bc_free(bc);
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
  ASSERT(mesh->num_cells + cond->bcell_map->size == N);
  for (int c = 0; c < mesh->num_cells; ++c)
    B[c] = mesh->cell_volumes[c] * cond->f[c];

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
  polymec_free(cond->f);
  polymec_free(cond->lambda);
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

krylov_solver_t* gmres_conduction_solver_new(mesh_t* mesh, int max_krylov_dim, int max_restarts)
{
  cond_t* cond = polymec_malloc(sizeof(cond_t));
  cond->mesh = mesh;
  cond->lambda = polymec_malloc(sizeof(real_t) * mesh->num_cells);
  cond->f = polymec_malloc(sizeof(real_t) * mesh->num_cells);
  for (int i = 0; i < mesh->num_cells; ++i)
  {
    cond->lambda[i] = 1.0;
    cond->f[i] = 0.0;
  }
  cond->bcs = string_ptr_unordered_map_new();
  cond->bcell_map = generate_boundary_cell_map(mesh);
  int N = mesh->num_cells + cond->bcell_map->size;
  krylov_solver_vtable vtable = {.ax = cond_ax, .b = cond_b, .dtor = cond_dtor};
  return gmres_krylov_solver_new(mesh->comm, cond, vtable, N, max_krylov_dim, max_restarts);
}

krylov_solver_t* bicgstab_conduction_solver_new(mesh_t* mesh, int max_krylov_dim)
{
  cond_t* cond = polymec_malloc(sizeof(cond_t));
  cond->mesh = mesh;
  cond->lambda = polymec_malloc(sizeof(real_t) * mesh->num_cells);
  cond->f = polymec_malloc(sizeof(real_t) * mesh->num_cells);
  for (int i = 0; i < mesh->num_cells; ++i)
  {
    cond->lambda[i] = 1.0;
    cond->f[i] = 0.0;
  }
  cond->bcs = string_ptr_unordered_map_new();
  cond->bcell_map = generate_boundary_cell_map(mesh);
  int N = mesh->num_cells + cond->bcell_map->size;
  krylov_solver_vtable vtable = {.ax = cond_ax, .b = cond_b, .dtor = cond_dtor};
  return bicgstab_krylov_solver_new(mesh->comm, cond, vtable, N, max_krylov_dim);
}

krylov_solver_t* tfqmr_conduction_solver_new(mesh_t* mesh, int max_krylov_dim)
{
  cond_t* cond = polymec_malloc(sizeof(cond_t));
  cond->mesh = mesh;
  cond->lambda = polymec_malloc(sizeof(real_t) * mesh->num_cells);
  cond->f = polymec_malloc(sizeof(real_t) * mesh->num_cells);
  for (int i = 0; i < mesh->num_cells; ++i)
  {
    cond->lambda[i] = 1.0;
    cond->f[i] = 0.0;
  }
  cond->bcs = string_ptr_unordered_map_new();
  cond->bcell_map = generate_boundary_cell_map(mesh);
  int N = mesh->num_cells + cond->bcell_map->size;
  krylov_solver_vtable vtable = {.ax = cond_ax, .b = cond_b, .dtor = cond_dtor};
  return tfqmr_krylov_solver_new(mesh->comm, cond, vtable, N, max_krylov_dim);
}

real_t* conduction_solver_lambda(krylov_solver_t* solver)
{
  cond_t* cond = krylov_solver_context(solver);
  return cond->lambda;
}

real_t* conduction_solver_f(krylov_solver_t* solver)
{
  cond_t* cond = krylov_solver_context(solver);
  return cond->f;
}

void conduction_solver_add_bc(krylov_solver_t* solver, const char* face_tag)
{
  cond_t* cond = krylov_solver_context(solver);

  // Retrieve the mesh's face tag.
  ASSERT(mesh_has_tag(cond->mesh->face_tags, face_tag));
  int num_faces;
  mesh_tag(cond->mesh->face_tags, face_tag, &num_faces);

  // Set up a boundary condition and stick it into our map of BCs.
  cond_bc_t* bc = cond_bc_new(num_faces);
  char* tag_name = string_dup(face_tag);
  string_ptr_unordered_map_insert_with_kv_dtor(cond->bcs, tag_name, bc, kv_dtor);
}

void conduction_solver_get_bc_arrays(krylov_solver_t* solver, const char* face_tag, 
                                     real_t** alpha, real_t** beta, real_t** gamma, int* num_faces)
{
  cond_t* cond = krylov_solver_context(solver);

  cond_bc_t** bc_ptr = (cond_bc_t**)string_ptr_unordered_map_get(cond->bcs, (char*)face_tag);
  if (bc_ptr != NULL)
  {
    cond_bc_t* bc = *bc_ptr;
    if (alpha != NULL)
      *alpha = bc->alpha;
    if (beta != NULL)
      *beta = bc->beta;
    if (gamma != NULL)
      *gamma = bc->gamma;
    if (num_faces != NULL)
      *num_faces = bc->num_faces;
  }
  else
  {
    if (alpha != NULL)
      *alpha = NULL;
    if (beta != NULL)
      *beta = NULL;
    if (gamma != NULL)
      *gamma = NULL;
    if (num_faces != NULL)
      *num_faces = 0;
  }
}

typedef struct
{
  mesh_t* mesh;
  real_t *lambda;
} cond_pc_t;

static void bjc_compute_diagonal(void* context, int block_size, 
                                 real_t alpha, real_t beta, real_t gamma, 
                                 real_t t, real_t* x, real_t* xdot, real_t* D)
{
  ASSERT(alpha == 0.0);
  ASSERT(beta == 1.0);
  ASSERT(gamma == 0.0);

  cond_pc_t* pc = context;
  mesh_t* mesh = pc->mesh;
  for (int cell = 0; cell < pc->mesh->num_cells; ++cell)
  {
    int pos = 0, face;

    // Evaluate the diagonal portion of div ( lambda o grad phi) using 
    // the discrete divergence theorem.
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
        grad = -x[cell] / point_distance(xn, xc);
        if (pc->lambda != NULL)
          lambda = 0.5 * (pc->lambda[cell] + pc->lambda[neighbor]);
      }
      else
      {
        point_t* xf = &mesh->face_centers[face];
        grad = x[cell] / (2.0 * point_distance(xf, xc));
        if (pc->lambda != NULL)
          lambda = pc->lambda[cell];
      }
      real_t area = mesh->face_areas[face];
      div_lambda_grad += area * lambda * grad;
    }
    D[cell] = div_lambda_grad;
  }
}

static void bjc_dtor(void* context)
{
  cond_pc_t* pc = context;
  polymec_free(pc->lambda);
  polymec_free(pc);
}

preconditioner_t* jacobi_conduction_pc_new(mesh_t* mesh)
{
  cond_pc_t* pc = polymec_malloc(sizeof(cond_pc_t));
  pc->mesh = mesh;
  pc->lambda = polymec_malloc(sizeof(real_t) * pc->mesh->num_cells);
  return block_jacobi_preconditioner_new(pc, bjc_compute_diagonal, bjc_dtor, pc->mesh->num_cells, 1);
}

real_t* jacobi_conduction_pc_lambda(preconditioner_t* pc)
{
  cond_pc_t* pcc = block_jacobi_preconditioner_context(pc);
  return pcc->lambda;
}

