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
#include "geometry/create_cvt.h"
#include "geometry/plane.h"
#include "geometry/cylinder.h"
#include "geometry/sphere.h"
#include "geometry/intersection.h"
#include "io/silo_io.h"
#include "io/vtk_plot_io.h"
#include "io/gnuplot_io.h"
#include "advect/advect_model.h"
#include "advect/advect_bc.h"
#include "advect/interpreter_register_advect_functions.h"
#include "advect/diffusion_op.h"

#ifdef __cplusplus
extern "C" {
#endif

// Linear solver context for the diffusion equation.
typedef struct
{
  KSP solver;
  Mat matrix;
  Vec solution;
  Vec rhs;
} advect_lin_sys_t;

static advect_lin_sys_t* advect_lin_sys_new(MPI_Comm comm, mesh_t* mesh, lin_op_t* op)
{
  advect_lin_sys_t* sys = malloc(sizeof(advect_lin_sys_t));
  MatCreate(comm, &sys->matrix);
  MatSetType(sys->matrix, MATSEQAIJ);
  MatSetOption(sys->matrix, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
//  MatSetOption(sys->matrix, MAT_NEW_NONZERO_LOCATION_ERR, PETSC_TRUE);
//  MatSetOption(sys->matrix, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);
  MatSetSizes(sys->matrix, mesh->num_cells, mesh->num_cells, PETSC_DETERMINE, PETSC_DETERMINE);
  VecCreate(comm, &sys->solution);
  VecSetType(sys->solution, VECSEQ);
  VecSetSizes(sys->solution, mesh->num_cells, PETSC_DECIDE);
  VecCreate(comm, &sys->rhs);
  VecSetType(sys->rhs, VECSEQ);
  VecSetSizes(sys->rhs, mesh->num_cells, PETSC_DECIDE);
  KSPCreate(comm, &sys->solver);

  // Pre-allocate matrix entries.
  PetscInt nz = 0, nnz[mesh->num_cells];
  for (int i = 0; i < mesh->num_cells; ++i)
    nnz[i] = lin_op_stencil_size(op, i);
  MatSeqAIJSetPreallocation(sys->matrix, nz, nnz);

  return sys;
}

static void advect_lin_sys_free(advect_lin_sys_t* sys)
{
  KSPDestroy(&sys->solver);
  MatDestroy(&sys->matrix);
  VecDestroy(&sys->solution);
  VecDestroy(&sys->rhs);
  free(sys);
}

// Advect model context structure.
typedef struct 
{
  mesh_t* mesh;             // Mesh.
  double* phi;              // Solution array.
  st_func_t* diffusivity;   // Diffusivity function.
  st_func_t* velocity;      // Velocity function.
  st_func_t* source;        // The (non-stiff, spatial) source term.
  st_func_t* solution;      // Analytic solution (if not NULL).
  st_func_t* initial_cond;  // The initial value of the solution.
  lin_op_t* D;              // Diffusion operator.

  str_ptr_unordered_map_t* bcs; // Boundary conditions.

  // Information for boundary cells.
  boundary_cell_map_t*  boundary_cells;

  // CFL safety factor.
  double CFL;             

  // This flag is true if the source function or one of the BCs is 
  // time-dependent, and false otherwise.
  bool is_time_dependent;

  // Diffusion equation linear system.
  advect_lin_sys_t* diff_system;

  // MPI communicator.
  MPI_Comm comm;            

} advect_t;

//------------------------------------------------------------------------
//                            Benchmarks
//------------------------------------------------------------------------

static model_t* create_advect(mesh_t* mesh,
                              st_func_t* velocity, 
                              st_func_t* diffusivity, 
                              st_func_t* source, 
                              st_func_t* initial_cond, 
                              str_ptr_unordered_map_t* bcs, 
                              st_func_t* solution,
                              options_t* options)
{
  // Create the model.
  model_t* model = advect_model_new(options);
  advect_t* a = (advect_t*)model_context(model);
  a->mesh = mesh;
  a->velocity = velocity;
  a->diffusivity = diffusivity;
  a->source = source;
  a->initial_cond = initial_cond;
  if (a->bcs != NULL)
    str_ptr_unordered_map_free(a->bcs);
  a->bcs = bcs;
  a->solution = solution;
  return model;
}

static void run_analytic_problem(model_t* model, 
                                 double t1, 
                                 double t2, 
                                 st_func_t* solution, 
                                 options_t* options,
                                 double* lp_norms)
{
  int max_steps = INT_MAX;
  char* max_steps_s = options_value(options, "max_steps");
  if (max_steps_s != NULL)
    max_steps = atoi(max_steps_s);

  // Run the thing.
  model_run(model, t1, t2, max_steps);

  // Calculate the Lp norm of the error and write it to Lp_norms.
  advect_t* a = model_context(model);
  double Linf = 0.0, L1 = 0.0, L2 = 0.0;
  for (int c = 0; c < a->mesh->num_cells; ++c)
  {
    double phi_sol;
    st_func_eval(solution, &a->mesh->cells[c].center, t2, &phi_sol);
    double V = a->mesh->cells[c].volume;
    double err = fabs(a->phi[c] - phi_sol);
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

static void advect_run_1d_flow(options_t* options, 
                               st_func_t* velocity, 
                               st_func_t* diffusivity, 
                               st_func_t* source, 
                               st_func_t* initial_cond, 
                               st_func_t* solution, 
                               int dim)
{
  // Boundary conditions.
  str_ptr_unordered_map_t* bcs = str_ptr_unordered_map_new();

  // Dirichlet on -x/+x.
  str_ptr_unordered_map_insert(bcs, "-x", advect_bc_new(1.0, 0.0, solution));
  str_ptr_unordered_map_insert(bcs, "+x", advect_bc_new(1.0, 0.0, solution));

  // Transverse faces - homogeneous Neumann BCs.
  double z = 0.0;
  st_func_t* zero = constant_st_func_new(1, &z);
  str_ptr_unordered_map_insert(bcs, "-y", advect_bc_new(0.0, 1.0, zero));
  str_ptr_unordered_map_insert(bcs, "+y", advect_bc_new(0.0, 1.0, zero));
  str_ptr_unordered_map_insert(bcs, "-z", advect_bc_new(0.0, 1.0, zero));
  str_ptr_unordered_map_insert(bcs, "+z", advect_bc_new(0.0, 1.0, zero));

  // Run times.
  double t1 = 0.0, t2 = 1.0;

  // Base resolution, number of runs.
  int N0;
  int num_runs;
  switch(dim)
  {
    case 1: 
      N0 = 32;
      num_runs = 4;
      break;
    case 2:
      N0 = 16;
      num_runs = 3;
      break;
    case 3:
      N0 = 8;
      num_runs = 3;
      break;
  }

  // Do a convergence study.
  double Lp_norms[num_runs][3];
  for (int iter = 0; iter < num_runs; ++iter)
  {
    int Nx = N0 * pow(2, iter), Ny = 1, Nz = 1;
    double dx = 1.0/Nx;
    bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
    if (dim == 1)
      bbox.y2 = bbox.z2 = dx;
    if (dim == 2)
    {
      Ny = Nx;
      bbox.z2 = dx;
    }
    if (dim == 3)
      Nz = Nx;
    mesh_t* mesh = create_cubic_lattice_mesh_with_bbox(Nx, Ny, Nz, &bbox);
    tag_cubic_lattice_mesh_faces(mesh, Nx, Ny, Nz, "-x", "+x", "-y", "+y", "-z", "+z");
    str_ptr_unordered_map_t* bcs_copy = str_ptr_unordered_map_copy(bcs);

    model_t* model = create_advect(mesh, velocity, diffusivity, source, 
                                   initial_cond, bcs_copy, solution, options);
    run_analytic_problem(model, t1, t2, solution, options, Lp_norms[iter]);
    model_free(model);

    // If we run in 1D or 2D, we need to adjust the norms.
    if (dim == 1)
    {
      Lp_norms[iter][1] *= Nx*Nx;
      Lp_norms[iter][2] *= Nx*Nx;
    }
    else if (dim == 2)
    {
      Lp_norms[iter][1] *= Nx;
      Lp_norms[iter][2] *= Nx;
    }
    log_urgent("iteration %d (Nx = %d): L1 = %g, L2 = %g, Linf = %g", iter, Nx, Lp_norms[iter][1], Lp_norms[iter][2], Lp_norms[iter][0]);
  }

  // Clean up.
  int pos = 0;
  char* key;
  void* value;
  while (str_ptr_unordered_map_next(bcs, &pos, &key, &value))
    advect_bc_free(value);
  str_ptr_unordered_map_free(bcs);
  zero = NULL;
}

static void run_stationary_flow_1d(options_t* options)
{
  double z = 0.0, o = 1.0;
  double v[] = {1.0, 0.0, 0.0};
  st_func_t* zero = constant_st_func_new(1, &z);
  st_func_t* one = constant_st_func_new(1, &o);
  st_func_t* v0 = constant_st_func_new(3, v);
  advect_run_1d_flow(options, v0, zero, zero, one, one, 1);
  zero = one = v0 = NULL;
}

static void stationary_blayer_1d_soln(void* ctx, point_t* x, double t, double* phi)
{
  double xx = x->x;
  double a = 1.0, d = 1.0;
  *phi = (exp(a/d) - exp(a*xx/d)) / (exp(a/d) - 1.0);
}

static void run_stationary_blayer_1d(options_t* options)
{
  double z = 0.0, o = 1.0;
  double v[] = {1.0, 0.0, 0.0};
  st_func_t* zero = constant_st_func_new(1, &z);
  st_func_t* one = constant_st_func_new(1, &o);
  st_func_t* v0 = constant_st_func_new(3, v);
  st_func_t* soln = st_func_from_func("blayer 1d", stationary_blayer_1d_soln,
                                      ST_INHOMOGENEOUS, ST_CONSTANT, 1);
  advect_run_1d_flow(options, v0, one, zero, soln, soln, 1);
  zero = one = v0 = soln = NULL;
}

static void square_wave_1d_soln(void* ctx, point_t* x, double t, double* phi)
{
  double width = 0.25;
  double a = 1.0;
  *phi = (fabs(x->x - a*t) < 0.5*width) ? 1.0 : 0.0;
}

static void run_square_wave_1d(options_t* options)
{
  double z = 0.0;
  double v[] = {1.0, 0.0, 0.0};
  st_func_t* zero = constant_st_func_new(1, &z);
  st_func_t* v0 = constant_st_func_new(3, v);
  st_func_t* soln = st_func_from_func("square wave 1d", square_wave_1d_soln,
                                      ST_INHOMOGENEOUS, ST_NONCONSTANT, 1);
  advect_run_1d_flow(options, v0, zero, zero, soln, soln, 1);
  zero = v0 = soln = NULL;
}

//------------------------------------------------------------------------
//                        Model implementation
//------------------------------------------------------------------------

// Compute half-state fluxes by upwinding (1st order).
static void compute_upwind_fluxes(mesh_t* mesh, st_func_t* velocity, st_func_t* source, boundary_cell_map_t* boundary_cells, double t, double dt, double* phi, double* fluxes)
{
  int num_faces = mesh->num_faces;
  for (int f = 0; f < num_faces; ++f)
  {
    face_t* face = &mesh->faces[f];
    if (face->cell2 == NULL) continue; // Boundary cell (handled below).

    // Compute the normal vector through this face, assuming that it is 
    // parallel to the displacement vector separating the centroids of the 
    // adjoining cells.
    vector_t n;
    point_displacement(&face->cell1->center, &face->cell2->center, &n);
    vector_normalize(&n);

    // Compute the normal velocity on this face at t + 0.5*dt.
    double v[3];
    st_func_eval(velocity, &face->center, t + 0.5*dt, v);
    double vn = v[0]*n.x + v[1]*n.y + v[2]*n.z;

    // If the normal velocity is positive, the velocity flows in the 
    // direction of cell 2 from cell 1, and cell 1 is upwind--otherwise cell2 
    // is upwind. In any case, we use the upwind value of the solution for 
    // the flux.
    int upwind_cell = (vn > 0.0) ? (face->cell1 - &mesh->cells[0]) 
                                 : (face->cell2 - &mesh->cells[0]);

    fluxes[f] = vn * phi[upwind_cell] * face->area;
// printf("%d,%d: vn = %g, F = %g\n", face->cell1 - &mesh->cells[0], face->cell2 - &mesh->cells[0], vn, fluxes[f]);

    // Cut off the fluxes below a certain threshold.
    if (fabs(fluxes[f]) < 1e-15) fluxes[f] = 0.0;
  }

  // Now go over the boundary cells and compute upwind fluxes.
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
    // ghost centroid. 
    for (int f = 0; f < cell_info->num_boundary_faces; ++f)
    {
      face_t* face = &mesh->faces[cell_info->boundary_faces[f]];
      int face_index = face - &mesh->faces[0];

      // Compute the normal vector through this face.
      vector_t n;
      point_displacement(&cell->center, &face->center, &n);
      vector_normalize(&n);

      // Compute the normal velocity on this face at t + 0.5*dt.
      double v[3];
      st_func_eval(velocity, &face->center, t + 0.5*dt, v);
      double vn = v[0]*n.x + v[1]*n.y + v[2]*n.z;

      // If the normal velocity is positive, the velocity flows in the 
      // direction of the ghost, and the interior cell 1 is upwind--otherwise 
      // the ghost cell is upwind. In any case, we use the upwind value of the 
      // solution for the flux.
      bool interior_is_upwind = (vn > 0.0);

      if (interior_is_upwind)
        fluxes[face_index] = vn * phi[bcell] * face->area;
      else
      {
        // The ghost cell is upwind, so we compute the solution there.

        // Retrieve the boundary condition for this face.
        advect_bc_t* bc = cell_info->bc_for_face[f];
        double alpha = bc->alpha, beta = bc->beta;

        // Compute L.
        double L = 2.0 * point_distance(&cell->center, &face->center);

        // Compute F at the face center.
        double F;
        st_func_eval(bc->F, &face->center, t, &F);

        // Compute the ghost value for the solution, and the resulting flux.
        double phi_g = (F + (beta/L - 0.5*alpha)) * phi[bcell] / (beta/L + 0.5*alpha);
// printf("%d: phi = %g, phi_g = %g\n", bcell, phi[bcell], phi_g);
        fluxes[face_index] = vn * phi_g * face->area;
      }
// printf("%d: vn = %g, F = %g\n", bcell, vn, fluxes[face_index]);
    }
  }
}

static void compute_half_step_fluxes(mesh_t* mesh, st_func_t* velocity, st_func_t* source, boundary_cell_map_t* boundary_cells, double t, double dt, double* phi, double* fluxes)
{
  // For now, we just use upwinding.
  compute_upwind_fluxes(mesh, velocity, source, boundary_cells, t, dt, phi, fluxes);
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
  }
  VecSetValues(b, mesh->num_cells, indices, values, INSERT_VALUES);
  VecAssemblyEnd(b);
}

static void compute_diffusive_deriv(mesh_t* mesh, 
                                    st_func_t* diffusivity,
                                    double* cell_sources,
                                    lin_op_t* D, 
                                    advect_lin_sys_t* diff_system,
                                    double t,
                                    double dt,
                                    double* diff_deriv)
{
#if 0
  // Set up the linear solver.
  KSPSetOperators(diff_system->solver, diff_system->matrix, diff_system->matrix, SAME_NONZERO_PATTERN);

  // Solve the linear system.
  KSPSolve(diff_system->solver, diff_system->rhs, diff_system->solution);

  // Access the solution values from the solution to our solution array.
  double* solution;
  VecGetArray(diff_system->solution, &solution);
  VecRestoreArray(diff_system->solution, &solution);
#endif

  // FIXME: For now, no diffusion.
  memset(diff_deriv, 0, sizeof(double)*mesh->num_cells);
}

static double advect_max_dt(void* context, double t, char* reason)
{
  advect_t* a = (advect_t*)context;

  // Find the minimum cell length / velocity ratio.
  double dt = FLT_MAX;
  for (int c = 0; c < a->mesh->num_cells; ++c)
  {
    // FIXME: For now, we just estimate the side of a cell from its 
    // FIXME: volume. This doesn't really work for nonuniform grids.
    double L = pow(a->mesh->cells[c].volume, 1.0/3.0);
    double V[3];
    st_func_eval(a->velocity, &a->mesh->cells[c].center, t, V);
    double Vmag = sqrt(V[0]*V[0] + V[1]*V[1] + V[2]*V[2]);
    if ((L / Vmag) < dt)
    {
      dt = a->CFL * L / Vmag;
      sprintf(reason, "CFL condition at cell %d.", c);
    }
  }

  return dt;
}

static void advect_advance(void* context, double t, double dt)
{
  advect_t* a = (advect_t*)context;

  // First, compute the half-step values of the fluxes on faces.
  double fluxes[a->mesh->num_faces];
  compute_half_step_fluxes(a->mesh, a->velocity, a->source, a->boundary_cells, t, dt, a->phi, fluxes);

  // Update the solution using the Divergence Theorem.
  int num_cells = a->mesh->num_cells;
  double phi_new[num_cells], diff_sources[num_cells];
  for (int c = 0; c < num_cells; ++c)
  {
    cell_t* cell = &a->mesh->cells[c];
    phi_new[c] = a->phi[c];
    for (int f = 0; f < cell->num_faces; ++f)
    {
      face_t* face = cell->faces[f];
      int face_index = face - &a->mesh->faces[0];

      // Make sure the sign of the flux is correct, since the fluxes have
      // been computed w.r.t. face->cell1.
      double sign = (face->cell1 == cell) ? 1.0 : -1.0;
      phi_new[c] -= sign * fluxes[face_index] * dt / cell->volume;
//printf("cell %d: flux %d (face %d) = %g\n", c, f, face_index, fluxes[face_index]);
    }

    // Add the non-stiff source term.
    double S;
    st_func_eval(a->source, &cell->center, t, &S);
    phi_new[c] += S * dt;

    // Compute the "source terms" that will be fed to the diffusion equation.
    double adv_deriv = (phi_new[c] - a->phi[c]) / dt;
    diff_sources[c] = S - adv_deriv;
  }

  // Do we have diffusivity? 
  if (a->diffusivity != NULL)
  {
    // Compute the diffusive derivative without splitting.
    double diff_deriv[num_cells];
    compute_diffusive_deriv(a->mesh, a->diffusivity, diff_sources, a->D, a->diff_system, t, dt, diff_deriv);

    // Update the solution.
    for (int c = 0; c < num_cells; ++c)
      phi_new[c] += diff_deriv[c] * dt;
  }

  // FIXME: Reactions go here.

  // Finally, update the model's solution vector.
  memcpy(a->phi, phi_new, sizeof(double)*num_cells);
}

static void advect_read_input(void* context, interpreter_t* interp, options_t* options)
{
  advect_t* a = (advect_t*)context;
  a->mesh = interpreter_get_mesh(interp, "mesh");
  if (a->mesh == NULL)
    polymec_error("advect: No mesh was specified.");
  a->velocity = interpreter_get_vector_function(interp, "velocity");
  if (a->velocity == NULL)
    polymec_error("advect: No velocity function was specified.");
  a->initial_cond = interpreter_get_scalar_function(interp, "initial_cond");
  if (a->initial_cond == NULL)
    polymec_error("advect: No initial condition (initial_cond) was specified.");
  a->source = interpreter_get_scalar_function(interp, "source");
  a->bcs = interpreter_get_table(interp, "bcs");
  if (a->bcs == NULL)
    polymec_error("poisson: No table of boundary conditions (bcs) was specified.");
  a->diffusivity = interpreter_get_scalar_function(interp, "diffusivity");

  // If CFL wasn't given on the command line, look for it here.
  if (options_value(options, "CFL") == NULL)
  {
    double CFL = interpreter_get_number(interp, "CFL");
    if (CFL != -FLT_MAX)
      a->CFL = CFL;
  }
}

static void advect_init(void* context, double t)
{
  advect_t* a = (advect_t*)context;
  ASSERT(a->mesh != NULL);
  ASSERT(a->initial_cond != NULL);
  ASSERT(a->velocity != NULL);
  ASSERT(st_func_num_comp(a->velocity) == 3);
  ASSERT(a->source != NULL);

  // Determine whether this model is time-dependent.
  a->is_time_dependent = !st_func_is_constant(a->source);
  int pos = 0;
  char* tag;
  advect_bc_t* bc;
  while (str_ptr_unordered_map_next(a->bcs, &pos, &tag, (void**)&bc))
  {
    // If any BC is time-dependent, the whole problem is.
    if (!st_func_is_constant(bc->F))
      a->is_time_dependent = true;
  }

  // Set up the diffusion operator.
  if ((a->diffusivity != NULL) && (a->D == NULL))
    a->D = diffusion_op_new(a->mesh, a->diffusivity);

  // If the model has been previously initialized, clean everything out.
  if (a->diff_system != NULL)
    advect_lin_sys_free(a->diff_system);
  if (a->phi != NULL)
    free(a->phi);

  if (a->boundary_cells != NULL)
    boundary_cell_map_free(a->boundary_cells);

  // Initialize the linear solver and friends.
  a->diff_system = advect_lin_sys_new(a->comm, a->mesh, a->D);

  // Gather information about boundary cells.
  a->boundary_cells = boundary_cell_map_from_mesh_and_bcs(a->mesh, a->bcs);

  // Initialize the solution.
  int num_cells = a->mesh->num_cells;
  a->phi = malloc(sizeof(double)*num_cells);
  for (int c = 0; c < num_cells; ++c)
    st_func_eval(a->initial_cond, &a->mesh->cells[c].center, t, &a->phi[c]);
}

static void advect_plot(void* context, io_interface_t* io, double t, int step)
{
  ASSERT(context != NULL);
  advect_t* a = (advect_t*)context;

  io_dataset_t* dataset = io_dataset_new("default");
  io_dataset_put_mesh(dataset, a->mesh);
  io_dataset_put_field(dataset, "phi", a->phi, 1, MESH_CELL, false);

  // Compute the cell-centered velocity and diffusivity and write them.
  double vel[a->mesh->num_cells], diff[a->mesh->num_cells];
  for (int c = 0; c < a->mesh->num_cells; ++c)
  {
    st_func_eval(a->velocity, &a->mesh->cells[c].center, t, &vel[c]);
    st_func_eval(a->diffusivity, &a->mesh->cells[c].center, t, &diff[c]);
  }
  io_dataset_put_field(dataset, "velocity", vel, 1, MESH_CELL, true);
  io_dataset_put_field(dataset, "diffusivity", diff, 1, MESH_CELL, true);

  // If we are given an analytic solution, write it and the solution error.
  if (a->solution != NULL)
  {
    int num_cells = a->mesh->num_cells;
    double soln[num_cells], error[num_cells];
    for (int c = 0; c < num_cells; ++c)
    {
      st_func_eval(a->solution, &a->mesh->cells[c].center, t, &soln[c]);
      error[c] = a->phi[c] - soln[c];
    }
    io_dataset_put_field(dataset, "solution", soln, 1, MESH_CELL, true);
    io_dataset_put_field(dataset, "error", error, 1, MESH_CELL, true);
  }

  io_append_dataset(io, dataset);
}

static void advect_save(void* context, io_interface_t* io, double t, int step)
{
  ASSERT(context != NULL);
  advect_t* a = (advect_t*)context;

  io_dataset_t* dataset = io_dataset_new("default");
  io_dataset_put_mesh(dataset, a->mesh);
  io_dataset_put_field(dataset, "phi", a->phi, 1, MESH_CELL, false);
  io_append_dataset(io, dataset);
}

static void advect_dtor(void* ctx)
{
  advect_t* a = (advect_t*)ctx;

  // Destroy BC table.
  str_ptr_unordered_map_free(a->bcs);

  if (a->mesh != NULL)
    mesh_free(a->mesh);
  if (a->diff_system != NULL)
    advect_lin_sys_free(a->diff_system);
  if (a->phi != NULL)
    free(a->phi);

  a->velocity = NULL;
  a->solution = NULL;
  a->initial_cond = NULL;
  a->diffusivity = NULL;
  a->source = NULL;

  boundary_cell_map_free(a->boundary_cells);
  free(a);
}

model_t* advect_model_new(options_t* options)
{
  model_vtable vtable = { .read_input = advect_read_input,
                          .init = advect_init,
                          .max_dt = advect_max_dt,
                          .advance = advect_advance,
                          .save = advect_save,
                          .plot = advect_plot,
                          .dtor = advect_dtor};
  advect_t* a = malloc(sizeof(advect_t));
  a->mesh = NULL;
  a->phi = NULL;
  a->diffusivity = NULL;
  a->velocity = NULL;
  double zero = 0.0;
  a->source = constant_st_func_new(1, &zero);
  a->solution = NULL;
  a->initial_cond = NULL;
  a->D = NULL;
  a->diff_system = NULL;
  a->bcs = str_ptr_unordered_map_new();
  a->boundary_cells = boundary_cell_map_new();
  a->comm = MPI_COMM_WORLD;
  model_t* model = model_new("advect", a, vtable, options);

  // Process options.
  {
    char* cfl_str = options_value(options, "CFL");
    if (cfl_str != NULL)
      a->CFL = atof(cfl_str);
    else
      a->CFL = 1.0;
    if ((a->CFL <= 0.0) || (a->CFL > 1.0))
      polymec_error("CFL should be between 0 and 1.");
  }

  // Register benchmarks.
  model_register_benchmark(model, "stationary_flow_1d", run_stationary_flow_1d, "Stational flow in 1D (v = 1).");
  model_register_benchmark(model, "stationary_blayer_1d", run_stationary_blayer_1d, "Stational flow with boundary layer in 1D.");
  model_register_benchmark(model, "square_wave_1d", run_square_wave_1d, "Square wave propagation in 1D.");

  // Set up saver/plotter.
  io_interface_t* saver = silo_io_new(MPI_COMM_SELF, 0, false);
  model_set_saver(model, saver);

  io_interface_t* plotter = NULL;
  char* which_plotter = options_value(options, "plotter");
  if (which_plotter != NULL)
  {
    if (!strcasecmp(which_plotter, "vtk"))
      plotter = vtk_plot_io_new(MPI_COMM_SELF, 0, false);
    else if (!strcasecmp(which_plotter, "silo"))
      plotter = silo_plot_io_new(MPI_COMM_SELF, 0, false);
    else if (!strcasecmp(which_plotter, "gnuplot"))
      plotter = gnuplot_io_new();
  }
  else
    plotter = vtk_plot_io_new(MPI_COMM_SELF, 0, false);
  if (plotter != NULL)
  {
    log_detail("Setting plotter to '%s'...", which_plotter);
    model_set_plotter(model, plotter);
  }

  return model;
}

#ifdef __cplusplus
}
#endif

