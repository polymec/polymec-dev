#include <string.h>
#include <stdlib.h>
#include "petscksp.h"
#include "petscmat.h"
#include "petscvec.h"
#include "core/unordered_map.h"
#include "core/least_squares.h"
#include "core/linear_algebra.h"
#include "core/constant_st_func.h"
#include "core/boundary_cell_map.h"
#include "geometry/cubic_lattice.h"
#include "geometry/create_cubic_lattice_mesh.h"
#include "geometry/interpreter_register_geometry_functions.h"
#include "io/silo_io.h"
#include "io/vtk_plot_io.h"
#include "poisson/poisson_model.h"
#include "poisson/poisson_bc.h"
#include "poisson/laplacian_op.h"
#include "poisson/interpreter_register_poisson_functions.h"
#include "poisson/register_poisson_benchmarks.h"

#ifdef __cplusplus
extern "C" {
#endif

// Poisson model context structure.
typedef struct 
{
  mesh_t* mesh;             // Mesh.
  st_func_t* rhs;           // Right-hand side function.
  double* phi;              // Solution array.
  lin_op_t* L;              // Laplacian operator.
  st_func_t* solution;      // Analytic solution (if non-NULL).

  str_ptr_unordered_map_t* bcs; // Boundary conditions.

  boundary_cell_map_t* boundary_cells; // Boundary cell info
  bool use_least_squares;   // Use least squares for boundary conditions?
  poly_ls_shape_t* shape;   // Least-squares shape functions.

  // This flag is true if the source function or one of the BCs is 
  // time-dependent, and false otherwise.
  bool is_time_dependent;

  KSP solver;               // Linear solver.
  Mat A;                    // Laplacian matrix.
  Vec x;                    // Solution vector.
  Vec b;                    // RHS vector. 

  bool initialized;         // Initialized flag.
  MPI_Comm comm;            // MPI communicator.

} poisson_t;

// Apply boundary conditions to a set of boundary cells using finite 
// differences.
static void apply_bcs_with_finite_differences(boundary_cell_map_t* boundary_cells,
                                              mesh_t* mesh,
                                              double t,
                                              Mat A,
                                              Vec b)
{
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
    // phi_g = (F + (beta/L - alpha/2)) * phi_i / (beta/L + alpha/2)
    //
    // where L is the distance between the interior centroid and the 
    // ghost centroid. These means that the only contributions to the 
    // linear system are on the diagonal of the matrix and to the 
    // right hand side vector.
    double Aii = 0.0, bi = 0;
    for (int f = 0; f < cell_info->num_boundary_faces; ++f)
    {
      // Retrieve the boundary condition for this face.
      poisson_bc_t* bc = cell_info->bc_for_face[f];

      face_t* face = &mesh->faces[cell_info->boundary_faces[f]];
      double alpha = bc->alpha, beta = bc->beta;

      // Compute L.
      double L = 2.0 * point_distance(&cell->center, &face->center);

      // Compute F at the face center.
      double F;
      st_func_eval(bc->F, &face->center, t, &F);

      // Add in the diagonal term (dphi/dn).
      Aii += ((beta/L - 0.5*alpha) / (beta/L + 0.5*alpha) - 1.0) * face->area / L;

      // Add in the right hand side contribution.
      bi -= (F / (beta/L + 0.5*alpha)) * face->area / L;
    }

    // Sum the values into the linear system.
    MatSetValues(A, 1, &bcell, 1, &bcell, &Aii, ADD_VALUES);
    VecSetValues(b, 1, &bcell, &bi, ADD_VALUES);
  }
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
  VecAssemblyEnd(b);
}

// Apply boundary conditions to a set of boundary cells using least 
// squares fits.
static void apply_bcs_with_least_squares(boundary_cell_map_t* boundary_cells,
                                         mesh_t* mesh,
                                         poly_ls_shape_t* shape,
                                         double t,
                                         Mat A,
                                         Vec b)
{
  // Go over the boundary cells, enforcing boundary conditions on each 
  // face.
  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  VecAssemblyBegin(b);
  int pos = 0, bcell;
  boundary_cell_t* cell_info;
  while (boundary_cell_map_next(boundary_cells, &pos, &bcell, &cell_info))
  {
    cell_t* cell = &mesh->cells[bcell];

    // Construct a polynomial least-squares fit for the cell and its 
    // neighbors (about its center), plus any boundary faces. That is,
    // construct a fit that satisfies the boundary conditions!

    // Number of points = number of neighbors + self + boundary faces
    int num_ghosts = cell_info->num_boundary_faces;
    int num_neighbors = cell_info->num_neighbor_cells;
    int num_points = num_neighbors + 1 + num_ghosts;
//printf("num_ghosts = %d, num_points = %d\n", num_ghosts, num_points);
    point_t points[num_points];
    point_copy(&points[0], &mesh->cells[bcell].center);
//printf("ipoint = %g %g %g\n", points[0].x, points[0].y, points[0].z);
    for (int n = 0; n < num_neighbors; ++n)
    {
      int neighbor = cell_info->neighbor_cells[n];
      point_copy(&points[n+1], &mesh->cells[neighbor].center);
    }
    int ghost_point_indices[num_ghosts];
    point_t constraint_points[num_ghosts];
    for (int n = 0; n < num_ghosts; ++n)
    {
      // We construct one ghost point per boundary face. This ghost point 
      // is the reflection of the boundary cell's centroid through the face.
      int bface = cell_info->boundary_faces[n];
      face_t* face = &mesh->faces[bface];
      int offset = 1 + num_neighbors;
      ghost_point_indices[n] = n+offset;
      points[n+offset].x = 2.0*face->center.x - cell->center.x;
      points[n+offset].y = 2.0*face->center.y - cell->center.y;
      points[n+offset].z = 2.0*face->center.z - cell->center.z;
//printf("gpoint = %g %g %g\n", points[n+offset].x, points[n+offset].y, points[n+offset].z);
      point_copy(&constraint_points[n], &face->center);
    }
    poly_ls_shape_set_domain(shape, &cell->center, points, num_points);

    // Now traverse the boundary faces of this cell and enforce the 
    // appropriate boundary condition on each. To do this, we compute 
    // the elements of an affine linear transformation that allows us 
    // to compute boundary values of the solution in terms of the 
    // interior values.
    vector_t face_normals[num_ghosts];
    double aff_matrix[num_ghosts*num_points], aff_vector[num_ghosts];
    {
      double a[num_ghosts], b[num_ghosts], c[num_ghosts], d[num_ghosts], e[num_ghosts];
      for (int f = 0; f < num_ghosts; ++f)
      {
        // Retrieve the boundary condition for the face.
        poisson_bc_t* bc = cell_info->bc_for_face[f];

        // Compute the face normal and the corresponding coefficients,
        // enforcing alpha * phi + beta * dphi/dn = F on the face.
        face_t* face = &mesh->faces[cell_info->boundary_faces[f]];
        vector_t* n = &face_normals[f];
        point_displacement(&cell->center, &face->center, n);
        vector_normalize(n);
        a[f] = bc->alpha;
        b[f] = bc->beta*n->x, c[f] = bc->beta*n->y, d[f] = bc->beta*n->z;

        // Compute F at the face center.
        st_func_eval(bc->F, &face->center, t, &e[f]);
      }

      // Now that we've gathered information about all the boundary
      // conditions, compute the affine transformation that maps the 
      // (unconstrained) values of the solution to the constrained values. 
      poly_ls_shape_compute_ghost_transform(shape, ghost_point_indices, num_ghosts, 
                                            constraint_points, a, b, c, d, e, aff_matrix, aff_vector);
//printf("%d: A_aff (num_ghosts = %d, np = %d) = ", bcell, num_ghosts, num_points);
//matrix_fprintf(aff_matrix, num_ghosts, num_points, stdout);
//printf("\n%d: b_aff = ", bcell);
//vector_fprintf(aff_vector, num_ghosts, stdout);
//printf("\n");
    }

    // Compute the flux through each boundary face and alter the 
    // linear system accordingly.
    int ij[num_neighbors+1];
    double N[num_points], Aij[num_neighbors+1];
    vector_t grad_N[num_points]; 
    for (int f = 0; f < num_ghosts; ++f)
    {
      // Compute the shape function values and gradients at the face center.
      face_t* face = &mesh->faces[cell_info->boundary_faces[f]];
      poly_ls_shape_compute_gradients(shape, &face->center, N, grad_N);
//printf("N = ");
//for (int i = 0; i < num_points; ++i)
//printf("%g ", N[i]);
//printf("\n");
//printf("grad N = ");
//for (int i = 0; i < num_points; ++i)
//printf("%g %g %g  ", grad_N[i].x, grad_N[i].y, grad_N[i].z);
//printf("\n");

      // Add the dphi/dn terms for face f to the matrix.
      vector_t* n = &face_normals[f];
      double bi = 0.0;

      // Diagonal term.
      ij[0] = bcell;
      // Compute the contribution to the flux from this cell.
//printf("For face %d (n = %g %g %g, x = %g %g %g):\n", f, n->x, n->y, n->z, face->center.x, face->center.y, face->center.z);
      Aij[0] = vector_dot(n, &grad_N[0]) * face->area; 
//printf("A[%d,%d] += %g * %g -> %g (%g)\n", bcell, ij[0], vector_dot(n, &grad_N[0]), face->area, Aij[0], N[0]);

      // Now compute the flux contributions from ghost points.
      for (int g = 0; g < num_ghosts; ++g)
      {
        double dNdn = vector_dot(n, &grad_N[num_neighbors+1+g]);
        Aij[0] += aff_matrix[num_ghosts*0+g] * dNdn * face->area;
//printf("A[%d,%d] += %g * %g * %g -> %g (%g)\n", bcell, ij[0], aff_matrix[g], dNdn, face->area, Aij[0], N[num_neighbors+1+g]);

        // Here we also add the affine contribution from this ghost
        // to the right-hand side vector.
        bi += -aff_vector[g] * dNdn * face->area;
      }

      // Compute the flux contributions from neighboring interior cells.
      for (int j = 0; j < num_neighbors; ++j)
      {
        ij[j+1] = cell_info->neighbor_cells[j];

        // Neighbor contribution.
        Aij[j+1] = vector_dot(n, &grad_N[j+1]) * face->area;

        // Ghost contributions.
        for (int g = 0; g < num_ghosts; ++g)
        {
          double dNdn = vector_dot(n, &grad_N[num_neighbors+1+g]);
          Aij[j+1] += aff_matrix[num_ghosts*(j+1)+g] * dNdn * face->area;
//printf("A[%d,%d] += %g * %g * %g = %g (%g)\n", bcell, ij[j+1], aff_matrix[num_ghosts*(j+1)+g],vector_dot(n, &grad_N[j+1]), face->area, Aij[j+1], N[j+1]);
        }
      }

      MatSetValues(A, 1, &bcell, num_neighbors+1, ij, Aij, ADD_VALUES);
      VecSetValues(b, 1, &bcell, &bi, ADD_VALUES);
    }
  }
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
  VecAssemblyEnd(b);
}

static void set_up_linear_system(mesh_t* mesh, lin_op_t* L, st_func_t* rhs, double t, Mat A, Vec b)
{
  PetscInt nnz[mesh->num_cells];
  for (int i = 0; i < mesh->num_cells; ++i)
    nnz[i] = lin_op_stencil_size(L, i);

  // Set matrix entries.
  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  for (int i = 0; i < mesh->num_cells; ++i)
  {
    int indices[nnz[i]];
    double values[nnz[i]];
    lin_op_compute_stencil(L, i, indices, values);
    for (int j = 0; j < nnz[i]; ++j)
      indices[j] += i;
    MatSetValues(A, 1, &i, nnz[i], indices, values, INSERT_VALUES);
  }
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

  // Compute the RHS vector.
  double values[mesh->num_cells];
  int indices[mesh->num_cells];
  VecAssemblyBegin(b);
  for (int c = 0; c < mesh->num_cells; ++c)
  {
    indices[c] = c;
    point_t xc = {.x = 0.0, .y = 0.0, .z = 0.0};
    st_func_eval(rhs, &xc, t, &values[c]);
    values[c] *= mesh->cells[c].volume;
  }
  VecSetValues(b, mesh->num_cells, indices, values, INSERT_VALUES);
  VecAssemblyEnd(b);
}

static void poisson_advance(void* context, double t, double dt)
{
  poisson_t* p = (poisson_t*)context;

  // If we're time-dependent, recompute the linear system here.
  if (!p->is_time_dependent)
  {
    set_up_linear_system(p->mesh, p->L, p->rhs, t+dt, p->A, p->b);
    if (p->use_least_squares)
      apply_bcs_with_least_squares(p->boundary_cells, p->mesh, p->shape, t+dt, p->A, p->b);
    else
      apply_bcs_with_finite_differences(p->boundary_cells, p->mesh, t+dt, p->A, p->b);
  }
//  MatView(p->A, PETSC_VIEWER_STDOUT_SELF);
//  VecView(p->b, PETSC_VIEWER_STDOUT_SELF);

  // Set up the linear solver.
  KSPSetOperators(p->solver, p->A, p->A, SAME_NONZERO_PATTERN);

  // Solve the linear system.
  KSPSolve(p->solver, p->b, p->x);

  // Copy the values from x to our solution array.
//  VecView(p->x, PETSC_VIEWER_STDOUT_SELF);
  double* x;
  VecGetArray(p->x, &x);
  memcpy(p->phi, x, sizeof(double)*p->mesh->num_cells);
  VecRestoreArray(p->x, &x);
}

static void poisson_read_input(void* context, interpreter_t* interp, options_t* options)
{
  poisson_t* p = (poisson_t*)context;
  p->mesh = interpreter_get_mesh(interp, "mesh");
  if (p->mesh == NULL)
    polymec_error("poisson: No mesh was specified.");
  p->rhs = interpreter_get_scalar_function(interp, "rhs");
  if (p->rhs == NULL)
    polymec_error("poisson: No right hand side (rhs) was specified.");
  p->bcs = interpreter_get_table(interp, "bcs");
  if (p->bcs == NULL)
    polymec_error("poisson: No table of boundary conditions (bcs) was specified.");

  // Set up everything else.
  p->L = laplacian_op_new(p->mesh);

  // Check the mesh to make sure it has all the tags mentioned in the 
  // bcs table.
  int pos = 0;
  char* tag;
  poisson_bc_t* bc;
  while (str_ptr_unordered_map_next(p->bcs, &pos, &tag, (void**)&bc))
  {
    // Retrieve the tag for this boundary condition.
    if (!mesh_has_tag(p->mesh->face_tags, tag))
      polymec_error("poisson: Face tag '%s' was not found in the mesh.", tag);
  }
}

static void poisson_init(void* context, double t)
{
  poisson_t* p = (poisson_t*)context;

  if (p->L == NULL)
    p->L = laplacian_op_new(p->mesh);

  // Determine whether this model is time-dependent.
  p->is_time_dependent = !st_func_is_constant(p->rhs);
  int pos = 0;
  char* tag;
  poisson_bc_t* bc;
  while (str_ptr_unordered_map_next(p->bcs, &pos, &tag, (void**)&bc))
  {
    // If any BC is time-dependent, the whole problem is.
    if (!st_func_is_constant(bc->F))
      p->is_time_dependent = true;
  }

  // If the model has been previously initialized, clean everything out.
  if (p->initialized)
  {
    KSPDestroy(&p->solver);
    MatDestroy(&p->A);
    VecDestroy(&p->x);
    VecDestroy(&p->b);
    free(p->phi);
    boundary_cell_map_free(p->boundary_cells);
    p->initialized = false;
  }

  // Initialize the linear solver and friends.
  MatCreate(p->comm, &p->A);
  MatSetType(p->A, MATSEQAIJ);
  MatSetOption(p->A, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
//  MatSetOption(p->A, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);
//  MatSetOption(p->A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);
  MatSetSizes(p->A, p->mesh->num_cells, p->mesh->num_cells, PETSC_DETERMINE, PETSC_DETERMINE);
  VecCreate(p->comm, &p->x);
  VecSetType(p->x, VECSEQ);
  VecSetSizes(p->x, p->mesh->num_cells, PETSC_DECIDE);
  VecCreate(p->comm, &p->b);
  VecSetType(p->b, VECSEQ);
  VecSetSizes(p->b, p->mesh->num_cells, PETSC_DECIDE);
  KSPCreate(p->comm, &p->solver);

  // Initialize the solution vector.
  p->phi = malloc(sizeof(double)*p->mesh->num_cells);

  // Pre-allocate matrix entries.
  PetscInt nz = 0, nnz[p->mesh->num_cells];
  for (int i = 0; i < p->mesh->num_cells; ++i)
    nnz[i] = lin_op_stencil_size(p->L, i);
  MatSeqAIJSetPreallocation(p->A, nz, nnz);

  // Gather information about boundary cells.
  p->boundary_cells = boundary_cell_map_from_mesh_and_bcs(p->mesh, p->bcs);

  // If we're independent of time, set up the linear system here.
  if (!p->is_time_dependent)
  {
    set_up_linear_system(p->mesh, p->L, p->rhs, t, p->A, p->b);
    if (p->use_least_squares)
      apply_bcs_with_least_squares(p->boundary_cells, p->mesh, p->shape, t, p->A, p->b);
    else
      apply_bcs_with_finite_differences(p->boundary_cells, p->mesh, t, p->A, p->b);
  }

  // We simply solve the problem for t = 0.
  poisson_advance((void*)p, t, 0.0);

  // We are now initialized.
  p->initialized = true;
}

static void poisson_plot(void* context, io_interface_t* io, double t, int step)
{
  ASSERT(context != NULL);
  poisson_t* p = (poisson_t*)context;

  io_dataset_t* dataset = io_dataset_new("default");
  io_dataset_put_mesh(dataset, p->mesh);
  io_dataset_put_field(dataset, "phi", p->phi, 1, MESH_CELL, false);

  // If we are given an analytic solution, write it and the solution error.
  if (p->solution != NULL)
  {
    double soln[p->mesh->num_cells], error[p->mesh->num_cells];
    for (int c = 0; c < p->mesh->num_cells; ++c)
    {
      st_func_eval(p->solution, &p->mesh->cells[c].center, t, &soln[c]);
      error[c] = p->phi[c] - soln[c];
    }
    io_dataset_put_field(dataset, "solution", soln, 1, MESH_CELL, true);
    io_dataset_put_field(dataset, "error", error, 1, MESH_CELL, true);
  }

#ifndef NDEBUG
  // If we're in debug mode, compute the laplacian of phi and dump that, too.
  // This ignores boundary conditions, so it's only useful as a diagnostic.
  double Lphi[p->mesh->num_cells];
  lin_op_apply(p->L, p->phi, Lphi);
  io_dataset_put_field(dataset, "L(phi)", Lphi, 1, MESH_CELL, true);
#endif

  io_append_dataset(io, dataset);
}

static void poisson_save(void* context, io_interface_t* io, double t, int step)
{
  ASSERT(context != NULL);
  poisson_t* p = (poisson_t*)context;

  io_dataset_t* dataset = io_dataset_new("default");
  io_dataset_put_mesh(dataset, p->mesh);
  io_dataset_put_field(dataset, "phi", p->phi, 1, MESH_CELL, false);
  io_append_dataset(io, dataset);
}

static void poisson_compute_error_norms(void* context, st_func_t* solution, double t, double* lp_norms)
{
  poisson_t* p = (poisson_t*)context;
  double Linf = 0.0, L1 = 0.0, L2 = 0.0;
  for (int c = 0; c < p->mesh->num_cells; ++c)
  {
    double phi_sol;
    st_func_eval(solution, &p->mesh->cells[c].center, t, &phi_sol);
    double V = p->mesh->cells[c].volume;
    double err = fabs(p->phi[c] - phi_sol);
//printf("i = %d, phi = %g, phi_s = %g, err = %g\n", c, a->phi[c], phi_sol, err);
    Linf = (Linf < err) ? err : Linf;
    L1 += err*V;
    L2 += err*err*V*V;
  }
  L2 = sqrt(L2);
  lp_norms[0] = Linf;
  lp_norms[1] = L1;
  lp_norms[2] = L2;
//  norm_t* lp_norm = fv2_lp_norm_new(a->mesh);
//  for (int p = 0; p <= 2; ++p)
//    lp_norms[p] = fv2_lp_norm_compute_error_from_solution(p, 
  // FIXME

  // Clean up.
//  lp_norm = NULL;
}

static void poisson_dtor(void* ctx)
{
  poisson_t* p = (poisson_t*)ctx;

  // Destroy BC table.
  str_ptr_unordered_map_free(p->bcs);

  if (p->mesh != NULL)
    mesh_free(p->mesh);
  if (p->initialized)
  {
    KSPDestroy(&p->solver);
    MatDestroy(&p->A);
    VecDestroy(&p->x);
    VecDestroy(&p->b);
    p->shape = NULL;
    free(p->phi);
  }

  boundary_cell_map_free(p->boundary_cells);
  free(p);
}

model_t* poisson_model_new(options_t* options)
{
  model_vtable vtable = { .read_input = poisson_read_input,
                          .init = poisson_init,
                          .advance = poisson_advance,
                          .save = poisson_save,
                          .plot = poisson_plot,
                          .compute_error_norms = poisson_compute_error_norms,
                          .dtor = poisson_dtor};
  poisson_t* context = malloc(sizeof(poisson_t));
  context->mesh = NULL;
  context->rhs = NULL;
  context->phi = NULL;
  context->L = NULL;
  context->bcs = str_ptr_unordered_map_new();
  context->solution = NULL;

  char* ls_opt = options_value(options, "use_least_squares");
  context->use_least_squares = false;
  if (ls_opt != NULL)
    context->use_least_squares = atoi(ls_opt);
  if (context->use_least_squares)
  {
    context->shape = poly_ls_shape_new(1, true);
    poly_ls_shape_set_simple_weighting_func(context->shape, 2, 1e-2);
  }
  else
    context->shape = NULL;
  context->boundary_cells = boundary_cell_map_new();
  context->initialized = false;
  context->comm = MPI_COMM_WORLD;
  model_t* model = model_new("poisson", context, vtable, options);

  // Set up an interpreter.
  interpreter_validation_t valid_inputs[] = {{"mesh", INTERPRETER_MESH},
                                             {"rhs", INTERPRETER_SCALAR_FUNCTION},
                                             {"bcs", INTERPRETER_TABLE},
                                             {"solution", INTERPRETER_SCALAR_FUNCTION},
                                             END_OF_VALID_INPUTS};
  model_enable_interpreter(model, valid_inputs);
  interpreter_register_geometry_functions(model_interpreter(model));
  interpreter_register_poisson_functions(model_interpreter(model));

  // Register benchmarks.
  register_poisson_benchmarks(model);

  // Set up saver/plotter.
  io_interface_t* saver = silo_io_new(MPI_COMM_SELF, 0, false);
  model_set_saver(model, saver);
  io_interface_t* plotter = vtk_plot_io_new(MPI_COMM_SELF, 0, false);
  model_set_plotter(model, plotter);

  return model;
}

model_t* create_poisson(mesh_t* mesh,
                        st_func_t* rhs,
                        str_ptr_unordered_map_t* bcs, 
                        st_func_t* solution,
                        options_t* options)
{
  model_t* model = poisson_model_new(options);
  poisson_t* pm = model_context(model);
  pm->mesh = mesh;
  pm->rhs = rhs;
  if (pm->bcs != NULL)
    str_ptr_unordered_map_free(pm->bcs);
  pm->bcs = bcs;
  pm->solution = solution;

  return model;
}

#ifdef __cplusplus
}
#endif

