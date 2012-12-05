#include "core/boundary_cell_map.h"
#include "advect/advect_diffusion_solver.h"
#include "advect/diffusion_op.h"
#include "advect/advect_bc.h"

#ifdef __cplusplus
extern "C" {
#endif

// Advective diffusion solver type
typedef struct
{
  st_func_t* diffusivity;
  lin_op_t* diff_op;
  st_func_t* source;
  mesh_t* mesh;
  boundary_cell_map_t* boundary_cells;
  double* advective_deriv;
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

  // Make sure the diffusion operator has the right diffusivity.
  diffusion_op_set_time(a->diff_op, t);

  // Now compute that thar matrix.
  PetscInt nnz[a->mesh->num_cells];
  for (int i = 0; i < a->mesh->num_cells; ++i)
    nnz[i] = lin_op_stencil_size(L, i);

  // Set matrix entries.
  MatAssemblyBegin(D, MAT_FLUSH_ASSEMBLY);
  for (int i = 0; i < a->mesh->num_cells; ++i)
  {
    int indices[nnz[i]];
    double values[nnz[i]];
    lin_op_compute_stencil(L, i, indices, values);
    for (int j = 0; j < nnz[i]; ++j)
      indices[j] += i;
    MatSetValues(D, 1, &i, nnz[i], indices, values, INSERT_VALUES);
  }
  MatAssemblyEnd(D, MAT_FLUSH_ASSEMBLY);

  // Now we have to set the entries for periodic boundaries.
  // We do this here so that the nonzero pattern of the matrix is properly
  // established.
  int pos = 0, bcell;
  boundary_cell_t* cell_info;
  MatAssemblyBegin(D, MAT_FINAL_ASSEMBLY);
  while (boundary_cell_map_next(a->boundary_cells, &pos, &bcell, &cell_info))
  {
    cell_t* cell = &a->mesh->cells[bcell];
    double Aij[2] = {0.0, 0.0};
    int ij[2];
    for (int f = 0; f < cell_info->num_boundary_faces; ++f)
    {
      // Retrieve the boundary condition for this face.
      advect_bc_t* bc = cell_info->bc_for_face[f];
      face_t* face = &a->mesh->faces[cell_info->boundary_faces[f]];
      double area = face->area;

      // Compute the diffusivity at the face center.
      double Di;
      st_func_eval(a->diffusivity, &face->center, t, &Di);

      if (bc == NULL) // periodic BC
      {
        ASSERT(cell_info->opp_faces[f] != NULL);
        face_t* other_face = cell_info->opp_faces[f];
        cell_t* other_cell = other_face->cell1;

        // Compute L.
        double L = point_distance(&cell->center, &face->center) + 
          point_distance(&other_face->center, &other_cell->center);

        // Add in the diagonal term (D * dphi/dn).
        ij[0] = bcell;
        Aij[0] -= Di * area / L;

        // Add in the off-diagonal term (D * dphi/dn).
        ij[1] = other_cell - &a->mesh->cells[0];
        Aij[1] = Di * area / L;
        MatSetValues(D, 1, &bcell, 2, ij, Aij, ADD_VALUES);
      }
    }
  }
  MatAssemblyEnd(D, MAT_FINAL_ASSEMBLY);
}

static void ad_compute_source_vector(void* context, Vec S, double t)
{
  ad_solver_t* a = (ad_solver_t*)context;
  mesh_t* mesh = a->mesh;

  // Compute the RHS vector.
  double values[mesh->num_cells];
  int indices[mesh->num_cells];
  VecAssemblyBegin(S);
  for (int c = 0; c < mesh->num_cells; ++c)
  {
    indices[c] = c;
    point_t xc = {.x = 0.0, .y = 0.0, .z = 0.0};
    st_func_eval(a->source, &xc, t, &values[c]);
    values[c] += a->advective_deriv[c]; // Add in advective derivative.
    values[c] *= mesh->cells[c].volume;
  }
  VecSetValues(S, mesh->num_cells, indices, values, INSERT_VALUES);
  VecAssemblyEnd(S);
}

static void ad_apply_bcs(void* context, Mat A, Vec b, double t)
{
  ad_solver_t* a = (ad_solver_t*)context;
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
    double Aii = 0.0, bi = 0;
    for (int f = 0; f < cell_info->num_boundary_faces; ++f)
    {
      // Retrieve the boundary condition for this face.
      advect_bc_t* bc = cell_info->bc_for_face[f];
      if (bc != NULL) // non-periodic BC
      {
        face_t* face = &mesh->faces[cell_info->boundary_faces[f]];
        double alpha = bc->alpha, beta = bc->beta;

        // Compute L.
        double L = 2.0 * point_distance(&cell->center, &face->center);

        // Compute F at the face center.
        double F;
        st_func_eval(bc->F, &face->center, t, &F);

        // Then the diagonal term is the part of the boundary condition
        // equation that multiplies phi_i. 
        double r = (beta/L - 0.5*alpha) / (beta/L + 0.5*alpha);
        Aii += 0.5 * alpha * (r + 1.0) + beta * (r - 1.0) / L;

        // Compute the right hand side contribution.
        double f = F / (0.5*alpha + beta/L);
        bi += (1.0 - (0.5*alpha + beta/L)*f) * F;
      }
    }

    // Insert the values into the linear system.
    MatSetValues(A, 1, &bcell, 1, &bcell, &Aii, ADD_VALUES);
    VecSetValues(b, 1, &bcell, &bi, ADD_VALUES);
  }
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
  VecAssemblyEnd(b);
}

static void ad_dtor(void* context)
{
  ad_solver_t* a = (ad_solver_t*)context;
  free(a->advective_deriv);
  free(a);
}

diffusion_solver_t* advect_diffusion_solver_new(st_func_t* diffusivity,
                                                st_func_t* source,
                                                mesh_t* mesh,
                                                boundary_cell_map_t* boundary_cells)
{
  ASSERT(diffusivity != NULL);
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
  a->diff_op = diffusion_op_new(mesh, diffusivity);
  a->diffusivity = diffusivity;
  a->source = source;
  a->mesh = mesh; // borrowed
  a->boundary_cells = boundary_cells; // borrowed
  a->advective_deriv = malloc(sizeof(double)*mesh->num_cells);
  memset(a->advective_deriv, 0, sizeof(double)*mesh->num_cells);
  diffusion_solver_t* solver = diffusion_solver_new("advective diffusion solver", a, vtable);
  return solver;
}

void advect_diffusion_solver_set_advective_deriv(diffusion_solver_t* solver,
                                                 double* advective_deriv)
{
  ad_solver_t* a = (ad_solver_t*)diffusion_solver_context(solver);
  memcpy(a->advective_deriv, advective_deriv, sizeof(double)*a->mesh->num_cells);
}

#ifdef __cplusplus
}
#endif

