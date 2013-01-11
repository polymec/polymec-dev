#include <string.h>

#include "cvode/cvode.h"
#include "cvode/cvode_spgmr.h"
#if USE_MPI
#include "nvector/nvector_parallel.h"
#else
#include "nvector/nvector_serial.h"
#endif

#include "core/unordered_map.h"
#include "core/unordered_set.h"
#include "core/boundary_cell_map.h"
#include "core/least_squares.h"
#include "io/silo_io.h"
#include "io/vtk_plot_io.h"
#include "io/gnuplot_io.h"
#include "cnav/cnav_implicit_model.h"
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
  void* cvode;                         // CVODE time integrator.
  int order;                           // Order of time integrator.
  N_Vector U;                          // Computed solution.
  st_func_t* source;                   // Source function
  st_func_t* initial_cond;             // Initial conditions.
  str_ptr_unordered_map_t* bcs;        // Boundary conditions.
  boundary_cell_map_t* boundary_cells; // Boundary cell data.
  double abs_tol, rel_tol;             // Absolute and relative tolerances.
  double CFL;                          // CFL safety factor.

  // Workspace stuff.
  double* ls_fits;                     // Least-squares fit coefficients.
} cnav_implicit_t;

static inline N_Vector N_Vector_new(MPI_Comm comm, int dim)
{
#ifdef USE_MPI
  // Add all the local dimensions to find tot_dim.
  int tot_dim;
  MPI_Allreduce(&dim, &tot_dim, 1, MPI_INT, MPI_SUM, comm);
  return N_VNew_Parallel(comm, dim, tot_dim);
#else
  return N_VNew_Serial(dim);
#endif
}

static inline void N_Vector_free(N_Vector v)
{
#ifdef USE_MPI
  N_VDestroy_Parallel(v);
#else
  N_VDestroy_Serial(v);
#endif
}

static inline double* N_Vector_data(N_Vector v)
{
#ifdef USE_MPI
  return NV_DATA_P(v);
#else
  return NV_DATA_S(v);
#endif
}

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
                            double* ns_fits,
                            vector_t* u_fits,
                            double* E_fits)
{
  int basis_dim = poly_ls_basis_size(order);
  // FIXME
}

static void compute_bulk_conserved_quantities(int order, 
                                              double* ls_fits, 
                                              int num_species, 
                                              double* masses,
                                              point_t* x,
                                              int cell_index, 
                                              double* rho, 
                                              vector_t* rhoU, 
                                              double* rhoE)
{
  *rho = 0.0;
  rhoU->x = rhoU->y = rhoU->z = 0.0;
  *rhoE = 0.0;
  int basis_dim = poly_ls_basis_size(order);
  double basis[basis_dim];
  compute_poly_ls_basis_vector(order, x, basis);

  double ns_fits[num_species*basis_dim], E_fits[basis_dim];
  vector_t u_fits[basis_dim];
  extract_ls_fits(order, ls_fits, cell_index, num_species, 
                  ns_fits, u_fits, E_fits);
  for (int i = 0; i < basis_dim; ++i)
  {
    for (int s = 0; s < num_species; ++s)
      *rho += ns_fits[num_species*i+s] * masses[s] * basis[i];
    rhoU->x += *rho * u_fits[i].x * basis[i];
    rhoU->y += *rho * u_fits[i].y * basis[i];
    rhoU->z += *rho * u_fits[i].z * basis[i];
    *rhoE += *rho * E_fits[i] * basis[i];
  }
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
                                        vector_t* F)
{
  int_unordered_set_t* computed_faces = int_unordered_set_new();

  int num_species = cnav_eos_num_species(eos);
  double masses[num_species];
  cnav_eos_get_masses(eos, masses);

  int num_comp = 4 + num_species;
  int p = 2; // FIXME: This is only 2nd-order accurate!
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

      // Compute the (bulk) conserved quantities for each of these cells 
      // at the face using the least-squares fits.
      double rho1, rho2, rhoE1, rhoE2;
      vector_t rhoU1, rhoU2;
      compute_bulk_conserved_quantities(p, ls_fits, num_species, masses, &face->center,
                                        c, &rho1, &rhoU1, &rhoE1);
      compute_bulk_conserved_quantities(p, ls_fits, num_species, masses, &face->center,
                                        ncell_index, &rho2, &rhoU2, &rhoE2);

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
                                     vector_t* F)
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
  double* u = N_Vector_data(U);
  cnav_implicit_t* cnav = (cnav_implicit_t*)context;
  mesh_t* mesh = cnav->mesh;
  cnav_eos_t* eos = cnav->eos;
  boundary_cell_map_t* boundary_cells = cnav->boundary_cells;
  int num_species = cnav_eos_num_species(eos);
  int num_comp = 4 + num_species;

  // Compute least-squares fit coefficients for each component of the 
  // solution.
  int order = cnav->order;
  double* ls_fits = cnav->ls_fits;
  compute_least_squares_fit(order, mesh, u, ls_fits);

  // Compute the nonlinear fluxes from hydrodynamic advection.
  vector_t hydro_fluxes[num_comp * mesh->num_faces];
  compute_hydrodynamic_fluxes(order, eos, mesh, boundary_cells, ls_fits, hydro_fluxes);

  // Compute the fluxes from diffusion processes.
  vector_t diff_fluxes[num_comp * mesh->num_faces];
  compute_diffusive_fluxes(order, eos, mesh, boundary_cells, ls_fits, diff_fluxes);

  // Now we use the discrete Divergence Theorem to update U_dot.
  // FIXME: This logic is only 2nd-order accurate!
  double* u_dot = N_Vector_data(U_dot);
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

static int compute_F_ale(double t, N_Vector U, N_Vector U_dot, void* context)
{
  double* u = N_Vector_data(U);
  cnav_implicit_t* cnav = (cnav_implicit_t*)context;
  mesh_t* mesh = cnav->mesh;
  cnav_eos_t* eos = cnav->eos;
  boundary_cell_map_t* boundary_cells = cnav->boundary_cells;
  int num_species = cnav_eos_num_species(eos);
  int num_comp = 7 + num_species;

  // Compute least-squares fit coefficients for each component of the 
  // solution.
  int order = cnav->order;
  double* ls_fits = cnav->ls_fits;
  compute_least_squares_fit(order, mesh, u, ls_fits);

  // Compute the nonlinear fluxes from hydrodynamic advection.
  vector_t hydro_fluxes[num_comp * mesh->num_faces];
  compute_hydrodynamic_fluxes(order, eos, mesh, boundary_cells, ls_fits, hydro_fluxes);

  // Compute the fluxes from diffusion processes.
  vector_t diff_fluxes[num_comp * mesh->num_faces];
  compute_diffusive_fluxes(order, eos, mesh, boundary_cells, ls_fits, diff_fluxes);

  // Now we use the discrete Divergence Theorem to update U_dot.
  // FIXME: This logic is only 2nd-order accurate!
  double* u_dot = N_Vector_data(U_dot);
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
  cnav_implicit_t* cnav = (cnav_implicit_t*)context;

  // Find the minimum cell length / velocity ratio.
  double dt = FLT_MAX;
  // FIXME
  return dt;
}

static void cnav_advance(void* context, double t, double dt)
{
  cnav_implicit_t* cnav = (cnav_implicit_t*)context;

  // Take a step.
  double t_actual;
  int status = CVode(cnav->cvode, t + dt, cnav->U, &t_actual, CV_NORMAL);
  ASSERT(status != CV_MEM_NULL);
  ASSERT(status != CV_NO_MALLOC);
  ASSERT(status != CV_ILL_INPUT);
  ASSERT(status != CV_LINIT_FAIL);
  ASSERT(status != CV_LSOLVE_FAIL);
  ASSERT(status != CV_RHSFUNC_FAIL);

  if ((status != CV_SUCCESS) && (status != CV_TSTOP_RETURN) && 
      (status != CV_ROOT_RETURN)) // && (status != CV_FIRST_RHSFUNC_FAIL))
  {
    switch(status)
    {
      case CV_TOO_CLOSE:
        polymec_error("dt (%g) is too small.", dt);
        break;
      case CV_TOO_MUCH_WORK:
        polymec_error("Advance took too many internal steps.");
        break;
      case CV_TOO_MUCH_ACC:
        polymec_error("The integrator could not achieve the desired accuracy.");
        break;
      case CV_ERR_FAILURE:
        polymec_error("Too many failures during one internal step.");
        break;
      case CV_CONV_FAILURE:
        polymec_error("Too many convergence failures.");
        break;
      case CV_REPTD_RHSFUNC_ERR:
        polymec_error("Too many errors in the right hand side function.");
        break;
      case CV_UNREC_RHSFUNC_ERR:
        polymec_error("Integrator failed to recover from a right hand side error.");
        break;
      case CV_RTFUNC_FAIL:
        polymec_error("Integrator could not find a root.");
        break;
      default:
        polymec_error("The integrator failed.");
        break;
    }
  }

  // If t_actual > t + dt, the output is still interpolated at the 
  // time t + dt, but the internal time is t_actual.
}

static void cnav_reconnect(void* context, mesh_t* new_mesh)
{
}

static void cnav_init(void* context, double t)
{
  cnav_implicit_t* cnav = (cnav_implicit_t*)context;

  // Make sure our solution vector is allocated.
  if (cnav->U != NULL)
  {
    N_Vector_free(cnav->U);
    CVodeFree(&cnav->cvode);
  }
  cnav->cvode = CVodeCreate(CV_BDF, CV_NEWTON);
  CVodeSetUserData(cnav->cvode, cnav);
  cnav->U = N_Vector_new(cnav->comm, cnav->mesh->num_cells);

  // Initialize the solution.
  double* U = N_Vector_data(cnav->U);
  for (int c = 0; c < cnav->mesh->num_cells; ++c)
    st_func_eval(cnav->initial_cond, &cnav->mesh->cells[c].center, t, &U[c]);
  CVodeInit(cnav->cvode, compute_F_eulerian, t, cnav->U);

  // We try GMRES with the preconditioner applied on the left, and 
  // using modified Gram-Schmidt orthogonalization.
  CVSpgmr(cnav->cvode, PREC_LEFT, 0);
  CVSpilsSetGSType(cnav->cvode, MODIFIED_GS);

  // Set the Jacobian-times-vector function.
  //CVSpilsSetJacTimesVecFn(cnav->cvode, JoV);

  // Set the preconditioner solve and setup functions.
  //CVSpilsSetPreconditioner(cvode_mem, Precond, PSolve);
}

static void cnav_save(void* context, io_interface_t* io, double t, int step)
{
  cnav_implicit_t* cnav = (cnav_implicit_t*)context;
  // FIXME
}

static void cnav_dtor(void* context)
{
  cnav_implicit_t* cnav = (cnav_implicit_t*)context;
  if (cnav->U != NULL)
    N_Vector_free(cnav->U);
  if (cnav->cvode != NULL)
    CVodeFree(&cnav->cvode);
  free(cnav);
}

model_t* cnav_implicit_model_new(int order,
                                 options_t* options)
{
  ASSERT(order >= 1);
  ASSERT(order <= 2);
//  ASSERT(order <= 5);

  model_vtable vtable = { .init = cnav_init,
                          .max_dt = cnav_max_dt,
                          .advance = cnav_advance,
                          .save = cnav_save,
                          .dtor = cnav_dtor};
  cnav_implicit_t* cnav = malloc(sizeof(cnav_implicit_t));
  model_t* model = model_new("Implicit compressible Navier-Stokes", cnav, vtable, options);

  // Initialize the bookkeeping structures.
  cnav->order = -1;
  cnav->comm = MPI_COMM_WORLD;
  cnav->cvode = NULL;
  cnav->U = NULL;

  // Set up the saver.
  io_interface_t* saver = silo_io_new(cnav->comm, 0, false);
  model_set_saver(model, saver);

  return model;
}

model_t* create_cnav_implicit(int order,
                              mesh_t* mesh,
                              cnav_eos_t* equation_of_state,
//                             reaction_network_t* reactions,
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
  model_t* model = cnav_implicit_model_new(order, options);
  cnav_implicit_t* cnav = (cnav_implicit_t*)model_context(model);
  cnav->order = order;
  cnav->mesh = mesh;
  cnav->source = source;
  cnav->initial_cond = initial_cond;
  // FIXME
  return model;
}

#ifdef __cplusplus
}
#endif

