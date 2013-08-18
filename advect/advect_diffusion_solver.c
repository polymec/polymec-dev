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

#include "core/boundary_cell_map.h"
#include "advect/advect_diffusion_solver.h"
#include "advect/diffusion_op.h"
#include "advect/advect_bc.h"

// Advective diffusion solver type
typedef struct
{
  st_func_t* diffusivity;
  lin_op_t* diff_op;
  st_func_t* source;
  mesh_t* mesh;
  boundary_cell_map_t* boundary_cells;
  double* advective_source;
} ad_solver_t;

static void ad_compute_diffusion_matrix(void* context, double_table_t* D, double t)
{
  ad_solver_t* a = (ad_solver_t*)context;
  lin_op_t* L = a->diff_op;
  ASSERT(L != NULL);

  // Make sure the diffusion operator has the right diffusivity.
  diffusion_op_set_time(a->diff_op, t);

  // Now compute that thar matrix.
  int nnz[a->mesh->num_cells];
  for (int i = 0; i < a->mesh->num_cells; ++i)
    nnz[i] = lin_op_stencil_size(L, i);

  // Set matrix entries.
  for (int i = 0; i < a->mesh->num_cells; ++i)
  {
    int indices[nnz[i]];
    double values[nnz[i]];
    lin_op_compute_stencil(L, i, indices, values);
    for (int j = 0; j < nnz[i]; ++j)
    {
      // Scale down by the cell volume.
      values[j] /= a->mesh->cells[i].volume;

      // Turn the indices into offsets.
      indices[j] += i;

//printf("D(%d, %d) = %g\n", i, indices[j], values[j]);
      double_table_insert(D, i, indices[j], values[j]);
    }
  }
}

static void ad_compute_source_vector(void* context, double* S, double t)
{
  ad_solver_t* a = (ad_solver_t*)context;
  mesh_t* mesh = a->mesh;
  ASSERT(a->source != NULL);

  // Compute the RHS vector.
  for (int c = 0; c < mesh->num_cells; ++c)
  {
    point_t xc = {.x = 0.0, .y = 0.0, .z = 0.0};
    st_func_eval(a->source, &xc, t, &S[c]);
    S[c] += a->advective_source[c]; // Add in advective source.
  }
}

static void ad_apply_bcs(void* context, double_table_t* A, double* b, double t)
{
  ad_solver_t* a = (ad_solver_t*)context;
  ASSERT(a->diffusivity != NULL);
  mesh_t* mesh = a->mesh;
  boundary_cell_map_t* boundary_cells = a->boundary_cells;

  // Go over the boundary cells, enforcing boundary conditions on each 
  // face.
  int pos = 0, bcell;
  boundary_cell_t* cell_info;
  while (boundary_cell_map_next(boundary_cells, &pos, &bcell, &cell_info))
  {
    cell_t* cell = &mesh->cells[bcell];
    double V = cell->volume;

    // We use ghost cells that are reflections of the boundary cells 
    // through the boundary faces. For each of these cells, the 
    // boundary condition alpha*phi + beta*dphi/dn = F assigns the 
    // ghost value
    //
    // phi_g = (F + (beta/L - alpha/2) * phi_i) / (beta/L + alpha/2)
    //
    // where L is the distance between the interior centroid and the 
    // ghost centroid. These means that the only contributions to the 
    // linear system are on the diagonal of the matrix and to the 
    // right hand side vector.
    for (int f = 0; f < cell_info->num_boundary_faces; ++f)
    {
      // Compute the diffusivity at the face center.
      face_t* face = &mesh->faces[cell_info->boundary_faces[f]];
      double Df;
      st_func_eval(a->diffusivity, &face->center, t, &Df);

      // Retrieve the boundary condition for this face.
      advect_bc_t* bc = cell_info->bc_for_face[f];

      if (bc != NULL) // non-periodic BC
      {
        double alpha = bc->alpha, beta = bc->beta;

        // Compute L.
        double L = 2.0 * point_distance(&cell->center, &face->center);

        // Compute F at the face center.
        double F;
        st_func_eval(bc->F, &face->center, t, &F);
//printf("F(%g, %g) = %g (%s)\n", face->center.x, t, F, st_func_name(bc->F));

        // Add in the diagonal term for (Df * dphi/dn)/V.
        double Aii = *double_table_get(A, bcell, bcell);
        Aii += Df * ((beta/L - 0.5*alpha) / (beta/L + 0.5*alpha) - 1.0) * face->area / (V * L);
        double_table_insert(A, bcell, bcell, Aii);

        // Add in the right hand side contribution.
        b[bcell] += -Df * (F / (beta/L + 0.5*alpha)) * face->area / (V * L);
//if (bc->alpha > 0.0)
//printf("A[%d,%d] = %g, b[%d] = %g\n", bcell, bcell, Aii, bcell, bi);
      }
      else  // periodic BC
      {
        double V = cell->volume;
        ASSERT(cell_info->opp_faces[f] != NULL);
        face_t* other_face = cell_info->opp_faces[f];
        cell_t* other_cell = other_face->cell1;

        // Compute L.
        double L = point_distance(&cell->center, &face->center) + 
          point_distance(&other_face->center, &other_cell->center);

        // Add in the diagonal term (D * dphi/dn).
        // Don't forget to scale down by the volume of the cell.
        double Aii = *double_table_get(A, bcell, bcell);
        Aii -= Df * face->area / (V * L);
        double_table_insert(A, bcell, bcell, Aii);

        // Add in the off-diagonal term (D * dphi/dn).
        int j = other_cell - &a->mesh->cells[0];
        double Aij = *double_table_get(A, bcell, j);
        Aij += Df * face->area / (V * L);
        double_table_insert(A, bcell, j, Aij);
      }
    }
  }
}

static void ad_dtor(void* context)
{
  ad_solver_t* a = (ad_solver_t*)context;
  free(a->advective_source);
  free(a);
}

diffusion_solver_t* advect_diffusion_solver_new(mesh_t* mesh,
                                                boundary_cell_map_t* boundary_cells,
                                                index_space_t* index_space)
{
  ASSERT(mesh != NULL);
  ASSERT(boundary_cells != NULL);

  diffusion_solver_vtable vtable = 
    { .compute_diffusion_matrix = ad_compute_diffusion_matrix,
      .compute_source_vector    = ad_compute_source_vector,
      .apply_bcs                = ad_apply_bcs,
      .dtor                     = ad_dtor
    };
  ad_solver_t* a = malloc(sizeof(ad_solver_t));
  a->diff_op = NULL;
  a->diffusivity = NULL;
  a->source = NULL;
  a->mesh = mesh; // borrowed
  a->boundary_cells = boundary_cells; // borrowed
  a->advective_source = malloc(sizeof(double)*mesh->num_cells);
  memset(a->advective_source, 0, sizeof(double)*mesh->num_cells);
  diffusion_solver_t* solver = diffusion_solver_new("advective diffusion solver", a, vtable, index_space);
  return solver;
}

void advect_diffusion_solver_set_diffusivity(diffusion_solver_t* solver, st_func_t* diffusivity)
{
  ASSERT(st_func_num_comp(diffusivity) == 1);
  ad_solver_t* a = (ad_solver_t*)diffusion_solver_context(solver);
  a->diffusivity = diffusivity;
  a->diff_op = diffusion_op_new(a->mesh, diffusivity);
}

void advect_diffusion_solver_set_source(diffusion_solver_t* solver, st_func_t* source)
{
  ASSERT(st_func_num_comp(source) == 1);
  ad_solver_t* a = (ad_solver_t*)diffusion_solver_context(solver);
  a->source = source;
}

void advect_diffusion_solver_set_advective_source(diffusion_solver_t* solver,
                                                  double* advective_source)
{
  ad_solver_t* a = (ad_solver_t*)diffusion_solver_context(solver);
  memcpy(a->advective_source, advective_source, sizeof(double)*a->mesh->num_cells);
}

