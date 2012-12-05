#include "core/boundary_cell_map.h"
#include "poisson/poisson_elliptic_solver.h"
#include "poisson/laplacian_op.h"
#include "poisson/poisson_bc.h"

#ifdef __cplusplus
extern "C" {
#endif

// Poisson elliptic solver type
typedef struct
{
  lin_op_t* op;
  st_func_t* source;
  mesh_t* mesh;
  boundary_cell_map_t* boundary_cells;
} p_solver_t;

static void p_create_matrix(void* context, Mat* mat)
{
  p_solver_t* p = (p_solver_t*)context;
  MatCreate(MPI_COMM_WORLD, mat);
  MatSetType(*mat, MATSEQAIJ);
  MatSetOption(*mat, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
//  MatSetOption(*mat, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);
//  MatSetOption(*mat, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);
  MatSetSizes(*mat, p->mesh->num_cells, p->mesh->num_cells, PETSC_DETERMINE, PETSC_DETERMINE);

  // Pre-allocate matrix entries.
  PetscInt nz = 0, nnz[p->mesh->num_cells];
  for (int i = 0; i < p->mesh->num_cells; ++i)
    nnz[i] = lin_op_stencil_size(p->op, i) + 1;
  MatSeqAIJSetPreallocation(*mat, nz, nnz);
}

static void p_create_vector(void* context, Vec* vec)
{
  p_solver_t* p = (p_solver_t*)context;
  VecCreate(MPI_COMM_WORLD, vec);
  VecSetType(*vec, VECSEQ);
  VecSetSizes(*vec, p->mesh->num_cells, PETSC_DECIDE);
}

static void p_create_ksp(void* context, KSP* ksp)
{
  KSPCreate(MPI_COMM_WORLD, ksp);
}

static void p_compute_operator_matrix(void* context, Mat D, double t)
{
  p_solver_t* p = (p_solver_t*)context;
  lin_op_t* L = p->op;

  // Now compute that thar matrix.
  PetscInt nnz[p->mesh->num_cells];
  for (int i = 0; i < p->mesh->num_cells; ++i)
    nnz[i] = lin_op_stencil_size(L, i);

  // Set matrix entries.
  MatAssemblyBegin(D, MAT_FINAL_ASSEMBLY);
  for (int i = 0; i < p->mesh->num_cells; ++i)
  {
    int indices[nnz[i]];
    double values[nnz[i]];
    lin_op_compute_stencil(L, i, indices, values);
    for (int j = 0; j < nnz[i]; ++j)
    {
      // Scale down by the cell volume.
      values[j] /= p->mesh->cells[i].volume;

      // Turn the indices into offsets.
      indices[j] += i;
    }
    MatSetValues(D, 1, &i, nnz[i], indices, values, INSERT_VALUES);
  }
  MatAssemblyEnd(D, MAT_FINAL_ASSEMBLY);
}

static void p_compute_source_vector(void* context, Vec S, double t)
{
  p_solver_t* p = (p_solver_t*)context;
  mesh_t* mesh = p->mesh;

  // Compute the RHS vector.
  double values[mesh->num_cells];
  int indices[mesh->num_cells];
  VecAssemblyBegin(S);
  for (int c = 0; c < mesh->num_cells; ++c)
  {
    indices[c] = c;
    point_t xc = {.x = 0.0, .y = 0.0, .z = 0.0};
    st_func_eval(p->source, &xc, t, &values[c]);
  }
  VecSetValues(S, mesh->num_cells, indices, values, INSERT_VALUES);
  VecAssemblyEnd(S);
}

static void p_apply_bcs(void* context, Mat A, Vec b, double t)
{
  p_solver_t* p = (p_solver_t*)context;
  mesh_t* mesh = p->mesh;
  boundary_cell_map_t* boundary_cells = p->boundary_cells;

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
      // Retrieve the boundary condition for this face.
      face_t* face = &mesh->faces[cell_info->boundary_faces[f]];
      poisson_bc_t* bc = cell_info->bc_for_face[f];

      if (bc != NULL) // non-periodic BC
      {
        double alpha = bc->alpha, beta = bc->beta;

        // Compute L.
        double L = 2.0 * point_distance(&cell->center, &face->center);

        // Compute F at the face center.
        double F;
        st_func_eval(bc->F, &face->center, t, &F);
//printf("F(%g, %g) = %g (%s)\n", face->center.x, t, F, st_func_name(bc->F));

        // Add in the diagonal term for (dphi/dn)/V.
        double Aii = ((beta/L - 0.5*alpha) / (beta/L + 0.5*alpha) - 1.0) * face->area / (V * L);

        // Add in the right hand side contribution.
        double bi = -(F / (beta/L + 0.5*alpha)) * face->area / (V * L);
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

        // Add in the diagonal term (dphi/dn).
        // Don't forget to scale down by the volume of the cell.
        int ij[2];
        double Aij[2];
        ij[0] = bcell;
        Aij[0] -= face->area / (V * L);

        // Add in the off-diagonal term (dphi/dn).
        ij[1] = other_cell - &p->mesh->cells[0];
        Aij[1] = face->area / (V * L);
        MatSetValues(A, 1, &bcell, 2, ij, Aij, ADD_VALUES);
      }
    }
  }
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
  VecAssemblyEnd(b);
//MatView(A, PETSC_VIEWER_STDOUT_SELF);
//VecView(b, PETSC_VIEWER_STDOUT_SELF);
}

static void p_dtor(void* context)
{
  p_solver_t* p = (p_solver_t*)context;
  free(p);
}

elliptic_solver_t* poisson_elliptic_solver_new(st_func_t* source,
                                               mesh_t* mesh,
                                               boundary_cell_map_t* boundary_cells)
{
  ASSERT(diffusivity != NULL);
  ASSERT(mesh != NULL);
  ASSERT(boundary_cells != NULL);

  elliptic_solver_vtable vtable = 
  { .create_matrix            = p_create_matrix,
    .create_vector            = p_create_vector,
    .create_ksp               = p_create_ksp,
    .compute_operator_matrix  = p_compute_operator_matrix,
    .compute_source_vector    = p_compute_source_vector,
    .apply_bcs                = p_apply_bcs,
    .dtor                     = p_dtor
  };
  p_solver_t* p = malloc(sizeof(p_solver_t));
  p->op = laplacian_op_new(mesh);
  p->source = source;
  p->mesh = mesh; // borrowed
  p->boundary_cells = boundary_cells; // borrowed
  elliptic_solver_t* solver = elliptic_solver_new("Poisson elliptic solver", p, vtable);
  return solver;
}

#ifdef __cplusplus
}
#endif

