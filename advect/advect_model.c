#include <string.h>
#include "core/unordered_map.h"
#include "core/least_squares.h"
#include "core/linear_algebra.h"
#include "core/constant_st_func.h"
#include "core/boundary_cell_map.h"
#include "io/silo_io.h"
#include "io/vtk_plot_io.h"
#include "io/gnuplot_io.h"
#include "advect/advect_model.h"
#include "advect/advect_bc.h"
#include "advect/advect_diffusion_solver.h"
#include "advect/interpreter_register_advect_functions.h"
#include "advect/register_advect_benchmarks.h"
#include "advect/slope_estimator.h"

#ifdef __cplusplus
extern "C" {
#endif

// Advect model context structure.
typedef struct 
{
  mesh_t* mesh;             // Mesh.
  double* phi;              // Solution array.
  st_func_t* diffusivity;   // Whole diffusivity function.
  st_func_t* velocity;      // Velocity function.
  st_func_t* source;        // The (non-stiff, spatial, whole) source term.
  st_func_t* solution;      // Analytic solution (if not NULL).
  st_func_t* initial_cond;  // The initial value of the solution.

  // Species-specific information
  int num_species; 
  st_func_t** species_diffusivities;
  st_func_t** species_sources;

  // Tiny optimizations.
  // stationary is set to true if V == 0
  // have_diffusivity is set to false if D == 0.
  bool stationary, have_diffusivity;

  str_ptr_unordered_map_t* bcs; // Boundary conditions.

  // Information for boundary cells.
  boundary_cell_map_t*  boundary_cells;

  // Slope estimator.
  slope_estimator_t* slope_estimator;

  // CFL safety factor.
  double CFL;             

  // Diffusion equation solver.
  diffusion_solver_t* diff_solver;

  // MPI communicator.
  MPI_Comm comm;            

} advect_t;

static double estimate_slope(slope_estimator_t* slope_est,
                             mesh_t* mesh,
                             double* phi,
                             int upwind_cell,
                             int downwind_cell,
                             boundary_cell_map_t* boundary_cells,
                             double t)
{
  // Gather the neighbors of the upwind cell and the corresponding
  // solution values.
  cell_t* up_cell = &mesh->cells[upwind_cell];
  double phi_up = phi[upwind_cell];
  cell_t* down_cell = &mesh->cells[downwind_cell];
  double phi_down = phi[downwind_cell];
  int num_other_points = up_cell->num_faces - 1;
  point_t other_points[num_other_points];
  double other_phi[num_other_points];
  int i = 0;
  for (int ff = 0; ff < up_cell->num_faces; ++ff)
  {
    face_t* face = up_cell->faces[ff];
    cell_t* other_cell = face_opp_cell(face, up_cell);
    if (other_cell == NULL)
    {
      // The upwind cell is a boundary cell. 
      boundary_cell_t* cell_info = *boundary_cell_map_get(boundary_cells, upwind_cell);

      // Which boundary face are we on?
      int bf = 0;
      while (&mesh->faces[cell_info->boundary_faces[bf]] != face) ++bf;
      ASSERT(bf < cell_info->num_boundary_faces);

      // What's the boundary condition on this face?
      advect_bc_t* bc = cell_info->bc_for_face[bf];
      point_t* face_center = &up_cell->faces[ff]->center;
      if (bc != NULL) 
      {
        // The other point is a ghost centroid, reflected over the face.
        other_points[i].x = 2.0*face_center->x - up_cell->center.x;
        other_points[i].y = 2.0*face_center->y - up_cell->center.y;
        other_points[i].z = 2.0*face_center->z - up_cell->center.z;

        // Compute the ghost value for the solution, and the resulting flux.
        double F;
        st_func_eval(bc->F, face_center, t, &F);
        double L = 2.0 * point_distance(&up_cell->center, &face->center);
        other_phi[i] = (F + (bc->beta/L - 0.5*bc->alpha)) * phi[upwind_cell] / (bc->beta/L + 0.5*bc->alpha);
      }
      else 
      {
        // The other point is a periodic point.
        ASSERT(cell_info->opp_faces[bf] != NULL);
        face_t* opp_face = cell_info->opp_faces[bf];
        cell_t* opp_cell = opp_face->cell1;
        other_points[i].x = face_center->x + (opp_cell->center.x - opp_face->center.x);
        other_points[i].y = face_center->y + (opp_cell->center.y - opp_face->center.y);
        other_points[i].z = face_center->z + (opp_cell->center.z - opp_face->center.z);
        int opp_cell_index = opp_cell - &mesh->cells[0];
        other_phi[i] = phi[opp_cell_index];
      }
      ++i;
    }
    else if (other_cell != down_cell)
    {
      point_copy(&other_points[i], &other_cell->center);
      other_phi[i] = phi[other_cell - &mesh->cells[0]];
      ++i;
    }
  }
  return slope_estimator_value(slope_est, &up_cell->center,
                               phi_up, &down_cell->center,
                               phi_down, other_points, 
                               other_phi, num_other_points);
}

static void compute_half_step_fluxes(mesh_t* mesh, 
                                     st_func_t* velocity, 
                                     st_func_t* source, 
                                     boundary_cell_map_t* boundary_cells, 
                                     double t, 
                                     double dt, 
                                     double* phi, 
                                     slope_estimator_t* slope_est,
                                     double* fluxes)
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
    double L = vector_mag(&n); // Grab the inter-cell spacing real quick.
    vector_normalize(&n);

    // Compute the normal velocity on this face at t + 0.5*dt.
    double v[3];
    st_func_eval(velocity, &face->center, t + 0.5*dt, v);
    double vn = v[0]*n.x + v[1]*n.y + v[2]*n.z;

    // If the normal velocity is positive, the velocity flows in the 
    // direction of cell 2 from cell 1, and cell 1 is upwind--otherwise cell2 
    // is upwind. In any case, we use the upwind value of the solution for 
    // the flux.
    int cell1 = face->cell1 - &mesh->cells[0],
        cell2 = face->cell2 - &mesh->cells[0];
    int upwind_cell = (vn > 0.0) ? cell1 : cell2;
    int downwind_cell = (vn > 0.0) ? cell2 : cell1;

    fluxes[f] = vn * phi[upwind_cell] * face->area;

    if (slope_est != NULL)
    {
      // Compute the piecewise-linear contribution to the integral of the 
      // advection equation. This is taken from Leveque's 1992 book, p. 185.
      double nu = vn * dt / L;

      // Estimate the slope through the upwind cell.
      double slope = estimate_slope(slope_est, mesh, phi, 
                                    upwind_cell, downwind_cell, boundary_cells, t);
//if (slope != 0.0)
//printf("slope = %g\n", slope);
//      fluxes[f] += 0.5 * vn * (1.0 - nu) * L * slope;
      fluxes[f] += 0.5 * vn * (SIGN(nu) - nu) * L * slope * face->area;
    }
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

      // Upwind and downwind values of phi.
      double phi_up, phi_down;

      // Retrieve the boundary condition for this face and use it to compute
      // the "ghost" value of phi.
      advect_bc_t* bc = cell_info->bc_for_face[f];
      double phi_g, L;
      if (bc != NULL) // Regular boundary condition.
      {
        double alpha = bc->alpha, beta = bc->beta;

        // Compute L, the centroid-to-centroid spacing.
        L = 2.0 * point_distance(&cell->center, &face->center);

        // Compute F at the face center.
        double F;
        st_func_eval(bc->F, &face->center, t, &F);

        // Compute the ghost value for the solution, and the resulting flux.
        phi_g = (F + (beta/L - 0.5*alpha)) * phi[bcell] / (beta/L + 0.5*alpha);
//printf("%d: phi = %g, phi_g = %g\n", bcell, phi[bcell], phi_g);
      }
      else // Periodic boundary condition.
      {
        ASSERT(cell_info->opp_faces[f] != NULL);
        face_t* opp_face = cell_info->opp_faces[f];
        cell_t* opp_cell = opp_face->cell1;
        int opp_cell_index = opp_cell - &mesh->cells[0];
        phi_g = phi[opp_cell_index];

        // Compute L, the centroid-to-centroid spacing.
        L = point_distance(&cell->center, &face->center) + 
            point_distance(&opp_cell->center, &opp_face->center);
      }

      if (vn > 0.0) 
      {
        phi_up = phi[bcell];
        phi_down = phi_g;
      }
      else
      {
        phi_up = phi_g;
        phi_down = phi[bcell];
      }

      // First-order flux computation.
      fluxes[face_index] = vn * phi_up * face->area;

      if (slope_est != NULL)
      {
        // Compute the piecewise-linear contribution to the integral of the 
        // advection equation. This is taken from Leveque's 1992 book, p. 185.
        double nu = vn * dt / L;
        double slope = 0.0; // FIXME minmod(phi_down - phi_up, phi_up - phi_down) / L;
        fluxes[f] += 0.5 * vn * (1.0 - nu) * L * slope;
      }
// printf("%d: vn = %g, F = %g\n", bcell, vn, fluxes[face_index]);
    }
  }
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
    if (a->stationary)
    {
      if (a->CFL*L < dt)
      {
        dt = a->CFL * L;
        sprintf(reason, "dt set to grid spacing at cell %d by diffusion.", c);
      }
    }
    else
    {
      double V[3];
      st_func_eval(a->velocity, &a->mesh->cells[c].center, t, V);
      double Vmag = sqrt(V[0]*V[0] + V[1]*V[1] + V[2]*V[2]);
      if ((a->CFL*L/Vmag) < dt)
      {
        dt = a->CFL * L / Vmag;
        sprintf(reason, "CFL condition at cell %d.", c);
      }
    }
  }

  return dt;
}

static void advect_advance(void* context, double t, double dt)
{
  advect_t* a = (advect_t*)context;
  int num_cells = a->mesh->num_cells;

  double phi_old[num_cells], phi_new[num_cells];
  for (int s = 0; s < a->num_species; ++s)
  {
    // Extract the species-specific solution.
    for (int c = 0; c < num_cells; ++c)
      phi_old[c] = a->phi[a->num_species*c+s];

    if (!a->stationary)
    {
      // First, compute the half-step values of the fluxes on faces.
      double fluxes[a->mesh->num_faces];
      compute_half_step_fluxes(a->mesh, a->velocity, a->source, a->boundary_cells, t, dt, phi_old, a->slope_estimator, fluxes);

      // Update the solution using the Divergence Theorem.
      double adv_source[num_cells];
      for (int c = 0; c < num_cells; ++c)
      {
        cell_t* cell = &a->mesh->cells[c];
        phi_new[c] = phi_old[c];
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
        adv_source[c] = -(phi_new[c] - phi_old[c]) / dt;
      }

      // Give the advective derivative to the diffusion solver.
      advect_diffusion_solver_set_advective_source(a->diff_solver, adv_source);
    }

    if (a->have_diffusivity)
    {
      // Set the species-specific diffusivity, source.
      advect_diffusion_solver_set_diffusivity(a->diff_solver, a->species_diffusivities[s]);
      advect_diffusion_solver_set_source(a->diff_solver, a->species_sources[s]);

      // Compute the diffusive derivative without splitting, using the 
      // 2nd-order L-stable TGA algorithm.
      diffusion_solver_tga(a->diff_solver, t, phi_old, t+dt, phi_new);
    }

    // Finally, update the model's solution vector.
    for (int c = 0; c < num_cells; ++c)
      a->phi[a->num_species*c+s] = phi_new[c];
  }
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
  ASSERT(a->diffusivity != NULL);
  ASSERT(st_func_num_comp(a->diffusivity) == st_func_num_comp(a->initial_cond));
  ASSERT(st_func_num_comp(a->source) == st_func_num_comp(a->initial_cond));
  ASSERT(a->velocity != NULL);
  ASSERT(st_func_num_comp(a->velocity) == 3);
  ASSERT(a->source != NULL);

  // Are we stationary (is V == 0)?
  a->stationary = true;
  if (st_func_is_constant(a->velocity) && 
      st_func_is_homogeneous(a->velocity))
  {
    double val[3];
    point_t x;
    st_func_eval(a->velocity, &x, 0.0, val);
    for (int i = 0; i < 3; ++i)
    {
      if (val[i] != 0.0)
      {
        a->stationary = false;
        break;
      }
    }
  }

  // Do we have nonzero diffusivity? 
  a->have_diffusivity = true;
  if (st_func_is_constant(a->diffusivity) && 
      st_func_is_homogeneous(a->diffusivity))
  {
    double val;
    point_t x;
    st_func_eval(a->diffusivity, &x, 0.0, &val);
    a->have_diffusivity = (val > 0.0);
  }

  // Extract species-specific information.
  a->num_species = st_func_num_comp(a->initial_cond);
  if (a->species_diffusivities != NULL)
  {
    free(a->species_diffusivities);
    free(a->species_sources);
  }
  a->species_diffusivities = malloc(sizeof(st_func_t*)*a->num_species);
  a->species_sources = malloc(sizeof(st_func_t*)*a->num_species);
  for (int s = 0; s < a->num_species; ++s)
  {
    a->species_diffusivities[s] = st_func_from_component(a->diffusivity, s);
    a->species_sources[s] = st_func_from_component(a->source, s);
  }

  // If the model has been previously initialized, clean everything out.
  if (a->diff_solver != NULL)
    diffusion_solver_free(a->diff_solver);
  if (a->phi != NULL)
    free(a->phi);

  if (a->boundary_cells != NULL)
    boundary_cell_map_free(a->boundary_cells);

  // Gather information about boundary cells.
  a->boundary_cells = boundary_cell_map_from_mesh_and_bcs(a->mesh, a->bcs);

  // Initialize the diffusion solver with an index space.
  index_space_t* is = index_space_new(a->comm, a->mesh->num_cells);
  a->diff_solver = advect_diffusion_solver_new(a->mesh, a->boundary_cells, is);

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
  for (int s = 0; s < a->num_species; ++s)
  {
    double phi[a->mesh->num_cells];
    for (int c = 0; c < a->mesh->num_cells; ++c)
      phi[c] = a->phi[a->num_species*c + s];
    char phi_name[128];
    sprintf(phi_name, "phi[%d]", s);
    io_dataset_put_field(dataset, phi_name, a->phi, 1, MESH_CELL, true);
  }

  // Compute the cell-centered velocity, diffusivities, and sources, 
  // and write them.
  double vel[a->mesh->num_cells];
  for (int c = 0; c < a->mesh->num_cells; ++c)
    st_func_eval(a->velocity, &a->mesh->cells[c].center, t, &vel[c]);
  io_dataset_put_field(dataset, "velocity", vel, 1, MESH_CELL, true);

  double diff[a->mesh->num_cells], source[a->mesh->num_cells];
  for (int s = 0; s < a->num_species; ++s)
  {
    for (int c = 0; c < a->mesh->num_cells; ++c)
    {
      st_func_eval(a->species_diffusivities[s], &a->mesh->cells[c].center, t, &diff[c]);
      st_func_eval(a->species_sources[s], &a->mesh->cells[c].center, t, &source[c]);
    }
    char diff_name[128];
    sprintf(diff_name, "diffusivity[%d]", s);
    io_dataset_put_field(dataset, diff_name, diff, 1, MESH_CELL, true);
    char source_name[128];
    sprintf(source_name, "source[%d]", s);
    io_dataset_put_field(dataset, source_name, source, 1, MESH_CELL, true);
  }

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

static void advect_compute_error_norms(void* context, st_func_t* solution, double t, double* lp_norms)
{
  advect_t* a = (advect_t*)context;
  double Linf = 0.0, L1 = 0.0, L2 = 0.0;
  for (int c = 0; c < a->mesh->num_cells; ++c)
  {
    double phi_sol;
    st_func_eval(solution, &a->mesh->cells[c].center, t, &phi_sol);
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

static void advect_dtor(void* ctx)
{
  advect_t* a = (advect_t*)ctx;

  // Destroy BC table.
  str_ptr_unordered_map_free(a->bcs);

  if (a->mesh != NULL)
    mesh_free(a->mesh);
  if (a->diff_solver != NULL)
    diffusion_solver_free(a->diff_solver);
  if (a->phi != NULL)
    free(a->phi);

  if (a->species_diffusivities != NULL)
  {
    free(a->species_diffusivities);
    free(a->species_sources);
  }

  a->velocity = NULL;
  a->solution = NULL;
  a->initial_cond = NULL;
  a->diffusivity = NULL;
  a->source = NULL;
  a->slope_estimator = NULL;

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
                          .compute_error_norms = advect_compute_error_norms,
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
  a->slope_estimator = slope_estimator_new(SLOPE_LIMITER_MINMOD);
  a->diff_solver = NULL;
  a->bcs = str_ptr_unordered_map_new();
  a->boundary_cells = boundary_cell_map_new();
  a->comm = MPI_COMM_WORLD;

  a->num_species = -1;
  a->species_diffusivities = NULL;
  a->species_sources = NULL;

  model_t* model = model_new("advect", a, vtable, options);

  // Process options.
  {
    char* cfl_str = options_value(options, "CFL");
    if (cfl_str != NULL)
      a->CFL = atof(cfl_str);
    else
      a->CFL = 0.9;
    if ((a->CFL <= 0.0) || (a->CFL > 1.0))
      polymec_error("CFL should be between 0 and 1.");

    char* limiter_str = options_value(options, "limiter");
    if (limiter_str != NULL)
    {
      if (!strcasecmp(limiter_str, "upwind"))
        a->slope_estimator = NULL;
      else
      {
        slope_limiter_t limiter;
        if (!strcasecmp(limiter_str, "lax-wendroff"))
          limiter = SLOPE_LIMITER_LAX_WENDROFF;
        else if (!strcasecmp(limiter_str, "beam-warming"))
          limiter = SLOPE_LIMITER_BEAM_WARMING;
        else if (!strcasecmp(limiter_str, "vanleer"))
          limiter = SLOPE_LIMITER_VAN_LEER;
        else if (!strcasecmp(limiter_str, "superbee"))
          limiter = SLOPE_LIMITER_SUPERBEE;
        else if (!strcasecmp(limiter_str, "minmod"))
          limiter = SLOPE_LIMITER_MINMOD;
        else
          polymec_error("Unknown limiter: %s", limiter_str);
        a->slope_estimator = slope_estimator_new(limiter);
      }
    }
  }

  // Register benchmarks.
  register_advect_benchmarks(model);

  // Set up saver/plotter.
  io_interface_t* saver = silo_io_new(a->comm, 0, false);
  model_set_saver(model, saver);

  io_interface_t* plotter = NULL;
  char* which_plotter = options_value(options, "plotter");
  if (which_plotter != NULL)
  {
    if (!strcasecmp(which_plotter, "vtk"))
      plotter = vtk_plot_io_new(a->comm, 0, false);
    else if (!strcasecmp(which_plotter, "silo"))
      plotter = silo_plot_io_new(a->comm, 0, false);
    else if (!strcasecmp(which_plotter, "gnuplot"))
      plotter = gnuplot_io_new();
  }
  else
    plotter = vtk_plot_io_new(a->comm, 0, false);
  if (plotter != NULL)
  {
    log_detail("Setting plotter to '%s'...", which_plotter);
    model_set_plotter(model, plotter);
  }

  return model;
}

model_t* create_advect(mesh_t* mesh,
                       st_func_t* velocity, 
                       st_func_t* diffusivity, 
                       st_func_t* source, 
                       st_func_t* initial_cond, 
                       str_ptr_unordered_map_t* bcs, 
                       st_func_t* solution,
                       options_t* options)
{
  ASSERT(mesh != NULL);
  ASSERT(velocity != NULL);
  ASSERT(st_func_num_comp(velocity) == 3);
  ASSERT(diffusivity != NULL);
  ASSERT(source != NULL);
  ASSERT(initial_cond != NULL);
  ASSERT(st_func_num_comp(diffusivity) == st_func_num_comp(initial_cond));
  ASSERT(st_func_num_comp(source) == st_func_num_comp(initial_cond));

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

#ifdef __cplusplus
}
#endif

