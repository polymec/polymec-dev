#include <string.h>

#include "core/unordered_map.h"
#include "core/unordered_set.h"
#include "core/boundary_cell_map.h"
#include "core/least_squares.h"
#include "io/silo_io.h"
#include "io/vtk_plot_io.h"
#include "io/gnuplot_io.h"
#include "cnav/cnav_model.h"
#include "cnav/cnav_bc.h"

#ifdef __cplusplus
extern "C" {
#endif

// This data structure represents all implicitly-integrated implementations
// of the compressible Navier-Stokes model.
typedef struct 
{
  MPI_Comm comm;
  mesh_t* mesh;                        // Computational mesh.
  cnav_eos_t* eos;                     // Equation of state.
  integrator_t* integrator;            // Time integrator.
  double* u;                           // Computed solution.
  st_func_t* source;                   // Source function
  st_func_t* initial_cond;             // Initial conditions.
  str_ptr_unordered_map_t* bcs;        // Boundary conditions.
  boundary_cell_map_t* boundary_cells; // Boundary cell data.
  double abs_tol, rel_tol;             // Absolute and relative tolerances.
  double CFL;                          // CFL safety factor.

  // Workspace stuff.
  double* ls_fits;                     // Least-squares fit coefficients.
} cnav_t;

static void compute_least_squares_fit(int order, 
                                      mesh_t* mesh, 
                                      double* U,
                                      double* ls_fit)
{
}

static void extract_ls_fits(int order,
                            double* ls_fits,
                            int cell_index,
                            int num_species,
                            double* c_fits,
                            vector_t* u_fits,
                            double* E_fits)
{
  int basis_dim = poly_ls_basis_size(order);
  // FIXME
}

static void solve_riemann_problem(double rho1, vector_t* u1, double p1,
                                  double rho2, vector_t* u2, double p2,
                                  double* riemann_rho, vector_t* riemann_u, double* riemann_p)
{
  // FIXME
  *riemann_rho = 0.5 * (rho1 + rho2);
  *riemann_p   = 0.5 * (p1 + p2);
  riemann_u->x = 0.5 * (u1->x + u2->x);
  riemann_u->y = 0.5 * (u1->y + u2->y);
  riemann_u->z = 0.5 * (u1->z + u2->z);
}

static void compute_hydrodynamic_fluxes(int order,
                                        cnav_eos_t* eos, 
                                        mesh_t* mesh, 
                                        boundary_cell_map_t* boundary_cells,
                                        double* ls_fits, 
                                        double* F)
{
  int_unordered_set_t* computed_faces = int_unordered_set_new();

  int num_species = cnav_eos_num_species(eos);
  int basis_dim = poly_ls_basis_size(order);
  double basis[basis_dim];

  int num_comp = 4 + num_species;
  int p = 1; // FIXME: This is only 2nd-order accurate!
  for (int c = 0; c < mesh->num_cells; ++c)
  {
    cell_t* cell = &mesh->cells[c];
    for (int f = 0; f < cell->num_faces; ++f)
    {
      face_t* face = cell->faces[f];
      int face_index = face - &mesh->faces[0];

      // Skip boundary faces.
      cell_t* ncell = face_opp_cell(face, cell);
      if (ncell == NULL) continue;
      int ncell_index = ncell - &mesh->cells[0];

      // Have we already the flux through this face?
      if (int_unordered_set_contains(computed_faces, face_index))
        continue;

      // Compute the face's normal vector.
      // FIXME: Wrong!
      vector_t normal;
      point_displacement(&cell->center, &ncell->center, &normal);
      vector_normalize(&normal);

      // Compute the (bulk) conserved quantities for each of these cells 
      // at the face using the least-squares fits.
      double rho1, rho2, rhoE1, rhoE2;
      vector_t rhoU1, rhoU2;
      {
        double c_fits[num_species*basis_dim], E_fits[basis_dim];
        vector_t u_fits[basis_dim];
        point_t x;

        rho1 = 0.0;
        rhoU1.x = rhoU1.y = rhoU1.z = 0.0;
        rhoE1 = 0.0;
        x.x = face->center.x - cell->center.x;
        x.y = face->center.y - cell->center.y;
        x.z = face->center.z - cell->center.z;
        compute_poly_ls_basis_vector(order, &x, basis);
        extract_ls_fits(order, ls_fits, c, num_species, 
                        c_fits, u_fits, E_fits);
        for (int i = 0; i < basis_dim; ++i)
        {
          for (int s = 0; s < num_species; ++s)
            rho1 += c_fits[num_species*i+s] * basis[i];
          rhoU1.x += rho1 * u_fits[i].x * basis[i];
          rhoU1.y += rho1 * u_fits[i].y * basis[i];
          rhoU1.z += rho1 * u_fits[i].z * basis[i];
          rhoE1 += rho1 * E_fits[i] * basis[i];
        }

        rho2 = 0.0;
        rhoU2.x = rhoU2.y = rhoU2.z = 0.0;
        rhoE2 = 0.0;
        x.x = face->center.x - ncell->center.x;
        x.y = face->center.y - ncell->center.y;
        x.z = face->center.z - ncell->center.z;
        compute_poly_ls_basis_vector(order, &x, basis);
        extract_ls_fits(order, ls_fits, ncell_index, num_species, 
                        c_fits, u_fits, E_fits);
        for (int i = 0; i < basis_dim; ++i)
        {
          for (int s = 0; s < num_species; ++s)
            rho2 += c_fits[num_species*i+s] * basis[i];
          rhoU2.x += rho2 * u_fits[i].x * basis[i];
          rhoU2.y += rho2 * u_fits[i].y * basis[i];
          rhoU2.z += rho2 * u_fits[i].z * basis[i];
          rhoE2 += rho2 * E_fits[i] * basis[i];
        }
      }

      // Transform the bulk conserved quantities into primitive variables
      // (rho, u, p).
      vector_t u1 = {.x = rhoU1.x/rho1, .y = rhoU1.y/rho1, .z = rhoU1.z/rho1};
      vector_t u2 = {.x = rhoU2.x/rho2, .y = rhoU2.y/rho2, .z = rhoU2.z/rho2};
      double eps1 = rhoE1/rho1 - 0.5*vector_dot(&u1, &u1);
      double T1 = cnav_eos_temperature(eos, rho1, eps1);
      double p1 = cnav_eos_pressure(eos, rho1, T1);
      double eps2 = rhoE2/rho2 - 0.5*vector_dot(&u2, &u2);
      double T2 = cnav_eos_temperature(eos, rho2, eps2);
      double p2 = cnav_eos_pressure(eos, rho2, T2);

      // Compute the solution to the Riemann problem.
      double riemann_rho, riemann_p;
      vector_t riemann_u;
      solve_riemann_problem(rho1, &u1, p1, rho2, &u2, p2, &riemann_rho, &riemann_u, &riemann_p);
      double riemann_un = vector_dot(&riemann_u, &normal);
      double riemann_eps = cnav_eos_specific_internal_energy(eos, riemann_rho, riemann_p);
      double riemann_rhoE = riemann_rho * (0.5 * vector_dot(&riemann_u, &riemann_u) + riemann_eps);

      // Construct the normal flux at the interface.
      for (int s = 0; s < num_species; ++s)
        F[num_comp*face_index+s] = riemann_rho * riemann_un;
      F[num_comp*face_index+num_species]   = riemann_rho * riemann_un * riemann_u.x + riemann_p;
      F[num_comp*face_index+num_species+1] = riemann_rho * riemann_un * riemann_u.y + riemann_p;
      F[num_comp*face_index+num_species+2] = riemann_rho * riemann_un * riemann_u.z + riemann_p;
      F[num_comp*face_index+num_species+3] = riemann_un * (riemann_rhoE + riemann_p);
    }
  }

  // Compute the fluxes on the boundary cells.
  int pos = 0, bcell;
  boundary_cell_t* cell_info;
  while (boundary_cell_map_next(boundary_cells, &pos, &bcell, &cell_info))
  {
    cell_t* cell = &mesh->cells[bcell];
    // FIXME
  }

  // Clean up.
  int_unordered_set_free(computed_faces);
}

static void compute_diffusive_fluxes(int order,
                                     cnav_eos_t* eos, 
                                     mesh_t* mesh, 
                                     boundary_cell_map_t* boundary_cells,
                                     double* ls_fit, 
                                     double* F)
{
  for (int c = 0; c < mesh->num_cells; ++c)
  {
    // FIXME
  }

  // Compute the derivatives on the boundary cells.
  int pos = 0, bcell;
  boundary_cell_t* cell_info;
  while (boundary_cell_map_next(boundary_cells, &pos, &bcell, &cell_info))
  {
    cell_t* cell = &mesh->cells[bcell];
  }
}

static int compute_F_eulerian(double t, N_Vector U, N_Vector U_dot, void* context)
{
  double* u = NULL; // FIXME: N_Vector_data(U);
  cnav_t* cnav = (cnav_t*)context;
  mesh_t* mesh = cnav->mesh;
  cnav_eos_t* eos = cnav->eos;
  boundary_cell_map_t* boundary_cells = cnav->boundary_cells;
  int num_species = cnav_eos_num_species(eos);
  int num_comp = 4 + num_species;

  // Compute least-squares fit coefficients for each component of the 
  // solution.
  int order = integrator_order(cnav->integrator);
  double* ls_fits = cnav->ls_fits;
  compute_least_squares_fit(order, mesh, u, ls_fits);

  // Compute the nonlinear fluxes from hydrodynamic advection.
  double hydro_fluxes[num_comp * mesh->num_faces];
  compute_hydrodynamic_fluxes(order, eos, mesh, boundary_cells, ls_fits, hydro_fluxes);

  // Compute the fluxes from diffusion processes.
  double diff_fluxes[num_comp * mesh->num_faces];
  compute_diffusive_fluxes(order, eos, mesh, boundary_cells, ls_fits, diff_fluxes);

  // Now we use the discrete Divergence Theorem to update U_dot.
  // FIXME: This logic is only 2nd-order accurate!
  double* u_dot = NULL; // FIXME: N_Vector_data(U_dot);
  for (int c = 0; c < mesh->num_cells; ++c)
  {
    cell_t* cell = &mesh->cells[c];
    double V = cell->volume;

    // Compute any source terms for this cell and start with them.
    st_func_eval(cnav->source, &cell->center, t, &u_dot[num_comp*c]);
    
    // Add in the fluxes.
    for (int f = 0; f < cell->num_faces; ++f)
    {
      face_t* face = cell->faces[f];
      int face_index = face - &mesh->faces[0];
      double A = face->area;
      vector_t normal;
      point_displacement(&cell->center, &face->center, &normal);
      vector_normalize(&normal);
      for (int i = 0; i < num_comp; ++i)
      {
        double Fn = hydro_fluxes[num_comp*face_index+i];
        double Gn = diff_fluxes[num_comp*face_index+i];
        u_dot[num_comp*c+i] -= (A/V) * (Fn + Gn);
      }
    }
  }

  return 0;
}

static int compute_F_ale(double t, N_Vector U, N_Vector U_dot, void* context)
{
  double* u = NULL; // FIXME: N_VGetArrayObject(U);
  cnav_t* cnav = (cnav_t*)context;
  mesh_t* mesh = cnav->mesh;
  cnav_eos_t* eos = cnav->eos;
  boundary_cell_map_t* boundary_cells = cnav->boundary_cells;
  int num_species = cnav_eos_num_species(eos);
  int num_comp = 7 + num_species;

  // Compute least-squares fit coefficients for each component of the 
  // solution.
  int order = integrator_order(cnav->integrator);
  double* ls_fits = cnav->ls_fits;
  compute_least_squares_fit(order, mesh, u, ls_fits);

  // Compute the nonlinear fluxes from hydrodynamic advection.
  double hydro_fluxes[num_comp * mesh->num_faces];
  compute_hydrodynamic_fluxes(order, eos, mesh, boundary_cells, ls_fits, hydro_fluxes);

  // Compute the fluxes from diffusion processes.
  double diff_fluxes[num_comp * mesh->num_faces];
  compute_diffusive_fluxes(order, eos, mesh, boundary_cells, ls_fits, diff_fluxes);

  // Now we use the discrete Divergence Theorem to update U_dot.
  // FIXME: This logic is only 2nd-order accurate!
  double* u_dot = NULL; // FIXME: N_Vector_data(U_dot);
  for (int c = 0; c < mesh->num_cells; ++c)
  {
    cell_t* cell = &mesh->cells[c];
    double V = cell->volume;

    // Compute any source terms for this cell and start with them.
    st_func_eval(cnav->source, &cell->center, t, &u_dot[num_comp*c]);
    
    // Add in the fluxes.
    for (int f = 0; f < cell->num_faces; ++f)
    {
      face_t* face = cell->faces[f];
      int face_index = face - &mesh->faces[0];
      double A = face->area;
      vector_t normal;
      point_displacement(&cell->center, &face->center, &normal);
      vector_normalize(&normal);
      for (int i = 0; i < num_comp; ++i)
      {
        double Fn = vector_dot(&hydro_fluxes[num_comp*face_index+i], &normal);
        double Gn = vector_dot(&diff_fluxes[num_comp*face_index+i], &normal);
        u_dot[num_comp*c+i] -= (A/V) * (Fn + Gn);
      }
    }
  }

  return 0;
}

static double cnav_max_dt(void* context, double t, char* reason)
{
  cnav_t* cnav = (cnav_t*)context;

  // Find the minimum cell length / velocity ratio.
  double dt = FLT_MAX;
  // FIXME
  return dt;
}

static void cnav_advance(void* context, double t, double dt)
{
  cnav_t* cnav = (cnav_t*)context;
  integrator_step(cnav->integrator, t, t + dt, cnav->u);
}

static void cnav_reconnect(void* context, mesh_t* new_mesh)
{
}

static void cnav_init(void* context, double t)
{
  // Initialize the solution.
  cnav_t* cnav = (cnav_t*)context;
  for (int c = 0; c < cnav->mesh->num_cells; ++c)
    st_func_eval(cnav->initial_cond, &cnav->mesh->cells[c].center, t, &cnav->u[c]);
}

static void cnav_save(void* context, io_interface_t* io, double t, int step)
{
  cnav_t* cnav = (cnav_t*)context;
  // FIXME
}

static void cnav_dtor(void* context)
{
  cnav_t* cnav = (cnav_t*)context;
  if (cnav->u != NULL)
    free(cnav->u);
  free(cnav);
}

model_t* cnav_model_new(options_t* options)
{
  ASSERT(order >= 1);
  ASSERT(order <= 2);
//  ASSERT(order <= 5);

  model_vtable vtable = { .init = cnav_init,
                          .max_dt = cnav_max_dt,
                          .advance = cnav_advance,
                          .save = cnav_save,
                          .dtor = cnav_dtor};
  cnav_t* cnav = malloc(sizeof(cnav_t));
  model_t* model = model_new("Compressible Navier-Stokes", cnav, vtable, options);

  // Initialize the bookkeeping structures.
  cnav->comm = MPI_COMM_WORLD;
  cnav->u = NULL;

  // Set up the saver.
  io_interface_t* saver = silo_io_new(cnav->comm, 0, false);
  model_set_saver(model, saver);

  return model;
}

model_t* create_cnav(integrator_type_t integrator,
                     int order,
                     mesh_t* mesh,
                     cnav_eos_t* equation_of_state,
//                   reaction_network_t* reactions,
                     st_func_t* source, 
                     st_func_t* initial_cond, 
                     str_ptr_unordered_map_t* bcs, 
                     st_func_t* solution,
                     options_t* options)
{
  ASSERT(order >= 1);
  ASSERT(order <= 5);
  ASSERT(mesh != NULL);
  ASSERT(equation_of_state != NULL);
  ASSERT(source != NULL);
  ASSERT(initial_cond != NULL);
  ASSERT(st_func_num_comp(source) == st_func_num_comp(initial_cond));

  // Create the model.
  model_t* model = cnav_model_new(options);
  cnav_t* cnav = (cnav_t*)model_context(model);
  cnav->integrator = NULL;
  cnav->mesh = mesh;
  cnav->source = source;
  cnav->initial_cond = initial_cond;
  // FIXME
  return model;
}

#ifdef __cplusplus
}
#endif

