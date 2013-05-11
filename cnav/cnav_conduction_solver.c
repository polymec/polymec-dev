#include "core/boundary_cell_map.h"
#include "cnav/cnav_conduction_solver.h"
#include "cnav/cnav_viscous_solver.h"
#include "cnav/cnav_bc.h"

// advective diffusion solver type
typedef struct
{
  st_func_t* diffusivity;
  lin_op_t* diff_op;
  st_func_t* source;
  mesh_t* mesh;
  boundary_cell_map_t* boundary_cells;
  double* advective_source;
} ad_solver_t;

static void ad_create_matrix(void* context, Mat* mat)
{
  ad_solver_t* a = (ad_solver_t*)context;
  MatCreate(MPI_COMM_WORLD, mat);
  MatSetType(*mat, MATSEQAIJ);
  MatSetOption(*mat, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
//  MatSetOption(*mat, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);
//  MatSetOption(*mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);
  MatSetSizes(*mat, a->mesh->num_cells, a->mesh->num_cells, PETSC_DETERMINE, PETSC_DETERMINE);

  // Pre-allocate matrix entries.
  PetscInt nz = 0, nnz[a->mesh->num_cells];
  for (int i = 0; i < a->mesh->num_cells; ++i)
    nnz[i] = lin_op_stencil_size(a->diff_op, i) + 1;
  MatSeqAIJSetPreallocation(*mat, nz, nnz);
}

static void ad_create_vector(void* context, Vec* vec)
{
  ad_solver_t* a = (ad_solver_t*)context;
  VecCreate(MPI_COMM_WORLD, vec);
  VecSetType(*vec, VECSEQ);
  VecSetSizes(*vec, a->mesh->num_cells, PETSC_DECIDE);
}

static void ad_create_ksp(void* context, KSP* ksp)
{
  KSPCreate(MPI_COMM_WORLD, ksp);
}

static void ad_compute_diffusion_matrix(void* context, Mat D, double t)
{
  ad_solver_t* a = (ad_solver_t*)context;
  lin_op_t* L = a->diff_op;
  ASSERT(L != NULL);

  // Make sure the diffusion operator has the right diffusivity.
  diffusion_op_set_time(a->diff_op, t);

  // Now compute that thar matrix.
  PetscInt nnz[a->mesh->num_cells];
  for (int i = 0; i < a->mesh->num_cells; ++i)
    nnz[i] = lin_op_stencil_size(L, i);

  // Set matrix entries.
  MatAssemblyBegin(D, MAT_FINAL_ASSEMBLY);
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
    }
    MatSetValues(D, 1, &i, nnz[i], indices, values, INSERT_VALUES);
  }
  MatAssemblyEnd(D, MAT_FINAL_ASSEMBLY);
}

static void ad_compute_source_vector(void* context, Vec S, double t)
{
  ad_solver_t* a = (ad_solver_t*)context;
  mesh_t* mesh = a->mesh;
  ASSERT(a->source != NULL);

  // Compute the RHS vector.
  double values[mesh->num_cells];
  int indices[mesh->num_cells];
  VecAssemblyBegin(S);
  for (int c = 0; c < mesh->num_cells; ++c)
  {
    indices[c] = c;
    point_t xc = {.x = 0.0, .y = 0.0, .z = 0.0};
    st_func_eval(a->source, &xc, t, &values[c]);
    values[c] += a->advective_source[c]; // Add in advective source.
  }
  VecSetValues(S, mesh->num_cells, indices, values, INSERT_VALUES);
  VecAssemblyEnd(S);
}

static void ad_apply_bcs(void* context, Mat A, Vec b, double t)
{
  ad_solver_t* a = (ad_solver_t*)context;
  ASSERT(a->diffusivity != NULL);
  mesh_t* mesh = a->mesh;
  boundary_cell_map_t* boundary_cells = a->boundary_cells;

  // Go over the boundary cells, enforcing boundary conditions on each 
  // face.
  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(b);
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
      cnav_bc_t* bc = cell_info->bc_for_face[f];

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
        double Aii = Df * ((beta/L - 0.5*alpha) / (beta/L + 0.5*alpha) - 1.0) * face->area / (V * L);

        // Add in the right hand side contribution.
        double bi = -Df * (F / (beta/L + 0.5*alpha)) * face->area / (V * L);
//if (bc->alpha > 0.0)
//printf("A[%d,%d] = %g, b[%d] = %g\n", bcell, bcell, Aii, bcell, bi);

        MatSetValues(A, 1, &bcell, 1, &bcell, &Aii, ADD_VALUES);
        VecSetValues(b, 1, &bcell, &bi, ADD_VALUES);
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
        int ij[2];
        double Aij[2];
        ij[0] = bcell;
        Aij[0] -= Df * face->area / (V * L);

        // Add in the off-diagonal term (D * dphi/dn).
        ij[1] = other_cell - &a->mesh->cells[0];
        Aij[1] = Df * face->area / (V * L);
        MatSetValues(A, 1, &bcell, 2, ij, Aij, ADD_VALUES);
      }
    }
  }
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
  VecAssemblyEnd(b);
//MatView(A, PETSC_VIEWER_STDOUT_SELF);
//VecView(b, PETSC_VIEWER_STDOUT_SELF);
}

static void ad_dtor(void* context)
{
  ad_solver_t* a = (ad_solver_t*)context;
  free(a->advective_source);
  free(a);
}

diffusion_solver_t* cnav_diffusion_solver_new(mesh_t* mesh,
                                                boundary_cell_map_t* boundary_cells)
{
  ASSERT(mesh != NULL);
  ASSERT(boundary_cells != NULL);

  diffusion_solver_vtable vtable = 
  { .create_matrix            = ad_create_matrix,
    .create_vector            = ad_create_vector,
    .create_ksp               = ad_create_ksp,
    .compute_diffusion_matrix = ad_compute_diffusion_matrix,
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
  diffusion_solver_t* solver = diffusion_solver_new("advective diffusion solver", a, vtable);
  return solver;
}

void cnav_diffusion_solver_set_diffusivity(diffusion_solver_t* solver, st_func_t* diffusivity)
{
  ASSERT(st_func_num_comp(diffusivity) == 1);
  ad_solver_t* a = (ad_solver_t*)diffusion_solver_context(solver);
  a->diffusivity = diffusivity;
  a->diff_op = diffusion_op_new(a->mesh, diffusivity);
}

void cnav_diffusion_solver_set_source(diffusion_solver_t* solver, st_func_t* source)
{
  ASSERT(st_func_num_comp(source) == 1);
  ad_solver_t* a = (ad_solver_t*)diffusion_solver_context(solver);
  a->source = source;
}

void cnav_diffusion_solver_set_advective_source(diffusion_solver_t* solver,
                                                  double* advective_source)
{
  ad_solver_t* a = (ad_solver_t*)diffusion_solver_context(solver);
  memcpy(a->advective_source, advective_source, sizeof(double)*a->mesh->num_cells);
}

