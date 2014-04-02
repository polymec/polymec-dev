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

#include "core/polymec.h"
#include "core/unordered_map.h"
#include "core/least_squares.h"
#include "core/linear_algebra.h"
#include "core/polynomial.h"
#include "core/write_silo.h"
#include "core/norms.h"
#include "geometry/interpreter_register_geometry_functions.h"
#include "integrators/nonlinear_integrator.h"
#include "integrators/lu_preconditioners.h"
#include "integrators/polyhedron_integrator.h"
#include "model/boundary_cell_map.h"
#include "model/polynomial_fit.h"
#include "poisson/poisson_model.h"
#include "poisson/poisson_bc.h"

// Equips the given interpreter with functions specific to the Poisson model.
extern void interpreter_register_poisson_functions(interpreter_t* interpreter);

// Registers the benchmarks for the Poisson model.
extern void register_poisson_benchmarks(model_t* model);

// Poisson model structure.
typedef struct 
{
  mesh_t* mesh;             
  adj_graph_t* graph;
  st_func_t* rhs;           // Right-hand side function.
  st_func_t* lambda;        // "Conduction" operator (symmetric tensor). 
  real_t* phi;              // Solution array.
  st_func_t* solution;      // Analytic solution (if non-NULL).
  st_func_t* initial_guess; // Initial guess (if non-NULL).

  string_ptr_unordered_map_t* bcs; // Boundary conditions.
  boundary_cell_map_t* boundary_cells; // Boundary cell -> BC mapping

  // Here we store the fluxes between cells/points, and keep track of 
  // which fluxes have already been computed. There are two fluxes for each 
  // face.
  real_t* face_fluxes;

  // Contributions from the source term (rhs).
  real_t* cell_sources;

  // Polynomial and polynomial fit.
  polynomial_fit_t* poly_fit;

  // Quadrature rules -- regular and "special" (to accommodate symmetry).
  polyhedron_integrator_t* poly_quad_rule;
  int_ptr_unordered_map_t* special_quad_rules;

  // Nonlinear solver that integrates Poisson's equation.
  nonlinear_integrator_t* solver;

  // Number of iterations for most recent nonlinear solve.
  int num_iterations;
} poisson_t;

static real_t poisson_advance(void* context, real_t max_dt, real_t t)
{
  poisson_t* p = context;
  nonlinear_integrator_solve(p->solver, t+max_dt, p->phi, &p->num_iterations);
  return max_dt;
}

static void poisson_read_input(void* context, interpreter_t* interp, options_t* options)
{
  poisson_t* p = context;
  p->mesh = interpreter_get_mesh(interp, "mesh");
  if (p->mesh == NULL)
    polymec_error("poisson: mesh is not specified.");
  p->rhs = interpreter_get_scalar_function(interp, "rhs");
  if (p->rhs == NULL)
    polymec_error("poisson: No valid right hand side (rhs) was specified.");
  p->lambda = interpreter_get_sym_tensor_function(interp, "lambda");
  if (p->lambda == NULL)
  {
    real_t ones[6] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    p->lambda = constant_st_func_new(6, ones);
  }
  p->bcs = interpreter_get_table(interp, "bcs");
  if (p->bcs == NULL)
    polymec_error("poisson: No table of boundary conditions (bcs) was specified.");

  // Check the mesh to make sure it has all the tags mentioned in the 
  // bcs table.
  int pos = 0;
  char* tag;
  poisson_bc_t* bc;
  while (string_ptr_unordered_map_next(p->bcs, &pos, &tag, (void**)&bc))
  {
    // Retrieve the tag for this boundary condition.
    if (!mesh_has_tag(p->mesh->face_tags, tag))
      polymec_error("poisson: Face tag '%s' was not found in the mesh.", tag);
  }
}

static int poisson_residual(void* context, real_t t, real_t* u, real_t* F)
{
  poisson_t* p = context;

  // Loop over all the cells and compute the fluxes for each one.
  for (int cell = 0; cell < p->mesh->num_cells; ++cell)
  {
    // Retrieve the quadrature rule for this cell.
    polyhedron_integrator_t* quad_rule = p->poly_quad_rule;
    polyhedron_integrator_t** special_rule = 
      (polyhedron_integrator_t**)int_ptr_unordered_map_get(p->special_quad_rules, cell);
    if (special_rule != NULL)
    {
      // This cell is handled by a special quadrature (for symmetry or 
      // special boundary conditions, for example).
      quad_rule = *special_rule;
    }
    else
      polyhedron_integrator_set_domain(quad_rule, p->mesh, cell);

    // Retrieve any boundary information for this cell.
    boundary_cell_t* bcell = NULL;
    {
      boundary_cell_t** bcell_ptr = (boundary_cell_t**)boundary_cell_map_get(p->boundary_cells, cell);
      if (bcell_ptr != NULL)
        bcell = *bcell_ptr;
    }

    // Find a least-squares fit for the solution in the vicinity of this cell.
    {
      point_t* x0 = &p->mesh->cell_centers[cell];
      polynomial_fit_reset(p->poly_fit, x0);

      // Self contribution.
      {
        real_t phi = u[cell];
        polynomial_fit_add_scatter_datum(p->poly_fit, 0, phi, x0);
      }

      int pos = 0, face;
      while (mesh_cell_next_face(p->mesh, cell, &pos, &face))
      {
        int local_face_index = pos - 1;
        int neighbor = mesh_face_opp_cell(p->mesh, face, cell);
        if (neighbor != -1)
        {
          // Contributions from neighboring cell.
          point_t* x = &p->mesh->cell_centers[neighbor];
          real_t phi = u[neighbor];
          polynomial_fit_add_scatter_datum(p->poly_fit, 0, phi, x);
        }
        else
        {
          // Get boundary condition information.
          ASSERT(bcell != NULL);
          poisson_bc_t* bc = bcell->bc_for_face[local_face_index];

          // Boundary quadrature contributions.
          int pos1 = 0;
          point_t xb;
          vector_t nb;
          real_t wb;
          while (polyhedron_integrator_next_surface_point(quad_rule,
                                                          local_face_index,
                                                          &pos1, 
                                                          &xb, &nb, &wb))
          {
            real_t alpha = bc->alpha, beta = bc->beta, gamma;
            st_func_eval(bc->F, &xb, t, &gamma);
//printf("at x = (%g, %g, %g): alpha = %g, beta = %g, gamma = %g\n", xb.x, xb.y, xb.z, alpha, beta, gamma);
            polynomial_fit_add_robin_bc(p->poly_fit, 0, alpha, beta, &nb, gamma, &xb);
          }
        }
      }

      // Solve the least squares system.
//if ((cell == 0) || (cell == p->mesh->num_cells - 1))
//{
//printf("%d: ", cell);
//polynomial_fit_fprintf(p->poly_fit, stdout);
//}
      polynomial_fit_compute(p->poly_fit);
//if ((cell == 0) || (cell == p->mesh->num_cells - 1))
//{
//polynomial_fit_fprintf(p->poly_fit, stdout);
//}
    }

    // Face fluxes.
    int fpos = 0, face;
    while (mesh_cell_next_face(p->mesh, cell, &fpos, &face))
    {
      real_t face_flux = 0.0;
      int local_face_index = fpos - 1;
      int which_cell = (cell == p->mesh->face_cells[2*face]) ? 0 : 1;

      // Go over the quadrature points in the face and compute the flux
      // at each point, accumulating the integral in face_flux.
      int qpos = 0;
      point_t xq;
      vector_t nq;
      real_t wq;
      while (polyhedron_integrator_next_surface_point(quad_rule, local_face_index, &qpos, &xq, &nq, &wq))
      {
        // Compute the gradient of the solution using the polynomial fit.
        vector_t grad_phi;
        polynomial_fit_eval_deriv(p->poly_fit, 1, 0, 0, &xq, &grad_phi.x);
        polynomial_fit_eval_deriv(p->poly_fit, 0, 1, 0, &xq, &grad_phi.y);
        polynomial_fit_eval_deriv(p->poly_fit, 0, 0, 1, &xq, &grad_phi.z);

        // Evaluate the conduction operator lambda.
        real_t lambda[6];
        st_func_eval(p->lambda, &xq, t, lambda);

        // Form the flux contribution.
        vector_t n_o_lambda;
        if (st_func_num_comp(p->lambda) == 1)
        {
          n_o_lambda.x = lambda[0]*nq.x;
          n_o_lambda.y = lambda[0]*nq.y;
          n_o_lambda.z = lambda[0]*nq.z;
        }
        else
        {
          n_o_lambda.x = nq.x * lambda[0] + nq.y * lambda[3] + nq.z * lambda[2];
          n_o_lambda.y = nq.x * lambda[1] + nq.y * lambda[2] + nq.z * lambda[4];
          n_o_lambda.z = nq.x * lambda[2] + nq.y * lambda[4] + nq.z * lambda[5];
        }
        face_flux += wq * vector_dot(&n_o_lambda, &grad_phi);
      }

      p->face_fluxes[2*face + which_cell] = face_flux;
//printf("flux(%d, %d) = %g\n", cell, local_face_index, face_flux);
    }

    // Source term (right hand side).
    {
      int qpos = 0;
      point_t xq;
      real_t wq;
      real_t source = 0;
      while (polyhedron_integrator_next_volume_point(quad_rule, &qpos, &xq, &wq))
      {
        // Compute the right hand side.
        real_t rhs;
        st_func_eval(p->rhs, &xq, t, &rhs);

        // Form the source contribution.
        source += wq * rhs;
      }
      p->cell_sources[cell] = source;
    }
  }

  // Loop over cells again, enforcing conservation and compute the residual.
  for (int cell = 0; cell < p->mesh->num_cells; ++cell)
  {
    // Enforce conservation by averaging all the face fluxes between 
    // neighboring cells.
    int fpos = 0, face;
    while (mesh_cell_next_face(p->mesh, cell, &fpos, &face))
    {
      int flux1_index = (cell == p->mesh->face_cells[2*face]) ? 2*face : 2*face+1;
      int flux2_index = (flux1_index == 2*face) ? 2*face+1 : 2*face;
      int other_cell = mesh_face_opp_cell(p->mesh, face, cell);
      if ((other_cell != -1) && (cell < other_cell))
      {
        // Retrieve the fluxes at face 1 and face 2.
        real_t flux1 = p->face_fluxes[flux1_index];
        real_t flux2 = p->face_fluxes[flux2_index];

        // Now average their magnitudes and make sure that they agree so 
        // that our fluxes are conservative.
        real_t avg_flux = 0.5 * (fabs(flux1) + fabs(flux2));
        int sign = (flux1 != 0.0) ? SIGN(flux1) : SIGN(flux2);
        p->face_fluxes[flux1_index] = sign * avg_flux;
        p->face_fluxes[flux2_index] = sign * avg_flux;
//printf("avg_flux(%d, %d) = %g\n", offset1, offset2, avg_flux);
      }
    }

    // Now compute the residual in this cell.
    F[cell] = p->cell_sources[cell];
    fpos = 0;
    while (mesh_cell_next_face(p->mesh, cell, &fpos, &face))
    {
      int which_cell = (cell == p->mesh->face_cells[2*face]) ? 0 : 1;
      F[cell] -= p->face_fluxes[2*face + which_cell];
    }
  }
//  log_debug("L2[F] = %g", l2_norm(F, p->mesh->num_cells));
//printf("F = [");
//for (int i = 0; i < p->mesh->num_cells; ++i)
//printf("%g ", F[i]);
//printf("]\n");

  return 0;
}

// This helper clears all the data structures that are tied to a discrete domain.
static void poisson_clear(poisson_t* p)
{
  if (p->graph != NULL)
    adj_graph_free(p->graph);

  if (p->face_fluxes != NULL)
  {
    free(p->face_fluxes);
    p->face_fluxes = NULL;
  }
  if (p->cell_sources != NULL)
  {
    free(p->cell_sources);
    p->cell_sources = NULL;
  }

  if (p->phi != NULL)
  {
    free(p->phi);
    p->phi = NULL;
  }

  if (p->boundary_cells != NULL)
  {
    boundary_cell_map_free(p->boundary_cells);
    p->boundary_cells = NULL;
  }

  if (p->poly_quad_rule != NULL)
  {
    polyhedron_integrator_free(p->poly_quad_rule);
    p->poly_quad_rule = NULL;
  }
  if (p->poly_fit != NULL)
  {
    polynomial_fit_free(p->poly_fit);
    p->poly_fit = NULL;
  }

  if (p->solver != NULL)
  {
    nonlinear_integrator_free(p->solver);
    p->solver = NULL;
  }
}

static void poisson_init(void* context, real_t t)
{
  poisson_t* p = context;

  // Clear temporal data structures.
  poisson_clear(p);

  int N = p->mesh->num_cells;
  p->graph = graph_from_mesh_cells(p->mesh);
  p->poly_quad_rule = midpoint_polyhedron_integrator_new();

  // Initialize the solution vector with the initial guess (or ones).
  p->phi = malloc(sizeof(real_t)*N);
  if (p->initial_guess != NULL)
  {
    // FIXME: Only 2nd-order accurate at the moment.
    for (int c = 0; c < N; ++c)
      st_func_eval(p->initial_guess, &p->mesh->cell_centers[c], t, &p->phi[c]);
  }
  else
  {
    for (int c = 0; c < N; ++c)
      p->phi[c] = 1.0;
  }

  // Gather information about boundary cells.
  p->boundary_cells = boundary_cell_map_from_mesh_and_bcs(p->mesh, p->bcs);

  // Initialize the nonlinear solver.
  nonlinear_integrator_vtable vtable = {.eval = poisson_residual};
  p->solver = bicgstab_nonlinear_integrator_new("Poisson", p, MPI_COMM_WORLD, N, vtable, LINE_SEARCH, 15);

  // Polynomial fit. Degree 1 for now.
  int poly_degree = 1;
  p->poly_fit = polynomial_fit_new(1, poly_degree);

  // For now, Use LU preconditioning with the same residual function.
  preconditioner_t* lu_precond = lu_preconditioner_new(p, poisson_residual, NULL, p->graph);
  nonlinear_integrator_set_preconditioner(p->solver, lu_precond);

  // Allocate storage for cell face fluxes.
  p->face_fluxes = malloc(sizeof(real_t)*2*p->mesh->num_faces);
  memset(p->face_fluxes, 0, sizeof(real_t)*2*p->mesh->num_faces);

  // Allocate storage for cell source contributions.
  p->cell_sources = malloc(sizeof(real_t)*N);
  memset(p->cell_sources, 0, sizeof(real_t)*N);

  // Now we simply solve the problem for the initial time.
  bool success = nonlinear_integrator_solve(p->solver, t, p->phi, &p->num_iterations);
  if (!success)
  {
    nonlinear_integrator_diagnostics_t diags;
    nonlinear_integrator_get_diagnostics(p->solver, &diags);
    nonlinear_integrator_diagnostics_fprintf(&diags, stdout);
    polymec_error("poisson: nonlinear solve failed.");
  }
printf("phi = [");
for (int i = 0; i < p->mesh->num_cells; ++i)
  printf("%g ", p->phi[i]);
printf("]\n");
}

static void poisson_plot(void* context, const char* prefix, const char* directory, real_t t, int step)
{
  ASSERT(context != NULL);
  poisson_t* p = context;

  string_ptr_unordered_map_t* cell_fields = string_ptr_unordered_map_new();
  string_ptr_unordered_map_insert(cell_fields, "phi", p->phi);

  // If we are given an analytic solution, write it and the solution error.
  // FIXME: Only 2nd-order accurate at the moment.
  if (p->solution != NULL)
  {
    real_t *soln = malloc(sizeof(real_t) * p->mesh->num_cells), 
           *error = malloc(sizeof(real_t) * p->mesh->num_cells);
    for (int c = 0; c < p->mesh->num_cells; ++c)
    {
      st_func_eval(p->solution, &p->mesh->cell_centers[c], t, &soln[c]);
      error[c] = p->phi[c] - soln[c];
    }
    string_ptr_unordered_map_insert_with_v_dtor(cell_fields, "solution", soln, DTOR(free));
    string_ptr_unordered_map_insert_with_v_dtor(cell_fields, "error", error, DTOR(free));
  }
  write_silo_mesh(MPI_COMM_SELF, prefix, directory, 0, 1, 0, p->mesh, cell_fields, 0.0);
}

static void poisson_save(void* context, const char* prefix, const char* directory, real_t t, int step)
{
  ASSERT(context != NULL);
  poisson_t* p = context;

  string_ptr_unordered_map_t* cell_fields = string_ptr_unordered_map_new();
  string_ptr_unordered_map_insert(cell_fields, "phi", p->phi);
  write_silo_mesh(MPI_COMM_SELF, prefix, directory, 0, 1, 0, p->mesh, cell_fields, 0.0);
}

static void poisson_compute_error_norms(void* context, st_func_t* solution, real_t t, real_t* lp_norms)
{
  poisson_t* p = context;
  real_t Linf = 0.0, L1 = 0.0, L2 = 0.0;
  for (int c = 0; c < p->mesh->num_cells; ++c)
  {
    real_t phi_sol;
    st_func_eval(solution, &p->mesh->cell_centers[c], t, &phi_sol);
    real_t V = p->mesh->cell_volumes[c];
    real_t err = fabs(p->phi[c] - phi_sol);
    //printf("i = %d, phi = %g, phi_s = %g, err = %g\n", c, a->phi[c], phi_sol, err);
    Linf = (Linf < err) ? err : Linf;
    L1 += err*V;
    L2 += err*err*V*V;
  }

  L2 = rsqrt(L2);
  lp_norms[0] = Linf;
  lp_norms[1] = L1;
  lp_norms[2] = L2;
}

static void poisson_dtor(void* context)
{
  poisson_t* p = context;

  // Clear temporal data structures.
  poisson_clear(p);

  // Destroy BC table.
  string_ptr_unordered_map_free(p->bcs);

  if (p->mesh != NULL)
    mesh_free(p->mesh);
  if (p->poly_fit != NULL)
    polynomial_fit_free(p->poly_fit);

  int_ptr_unordered_map_free(p->special_quad_rules);

  if (p->phi != NULL)
    free(p->phi);

  if (p->boundary_cells != NULL)
    boundary_cell_map_free(p->boundary_cells);

  if (p->solver != NULL)
    nonlinear_integrator_free(p->solver);

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
  poisson_t* p = malloc(sizeof(poisson_t));
  p->mesh = NULL;
  p->graph = NULL;
  p->rhs = NULL;
  p->lambda = NULL;
  p->phi = NULL;
  p->solution = NULL;
  p->initial_guess = NULL;
  p->bcs = string_ptr_unordered_map_new();
  p->boundary_cells = boundary_cell_map_new();
  p->face_fluxes = NULL;
  p->cell_sources = NULL;
  p->poly_fit = NULL;
  p->poly_quad_rule = NULL;
  p->special_quad_rules = int_ptr_unordered_map_new();
  p->solver = NULL;

  model_t* model = model_new("poisson", p, vtable, options);

  // Set up an interpreter.
  interpreter_validation_t valid_inputs[] = {{"mesh", INTERPRETER_MESH},
                                             {"points", INTERPRETER_POINT_LIST},
                                             {"lambda", INTERPRETER_SYM_TENSOR_FUNCTION},
                                             {"rhs", INTERPRETER_SCALAR_FUNCTION},
                                             {"bcs", INTERPRETER_TABLE},
                                             {"solution", INTERPRETER_SCALAR_FUNCTION},
                                             END_OF_VALID_INPUTS};
  model_enable_interpreter(model, valid_inputs);
  interpreter_register_geometry_functions(model_interpreter(model));
  interpreter_register_poisson_functions(model_interpreter(model));

  // Register benchmarks.
  register_poisson_benchmarks(model);

  return model;
}

void poisson_model_set_initial_guess(model_t* model,
                                     st_func_t* guess)
{
  poisson_t* p = model_context(model);
  ASSERT(st_func_num_comp(guess) == 1);
  p->initial_guess = guess;
}

model_t* create_poisson(mesh_t* mesh,
                        st_func_t* lambda,
                        st_func_t* rhs,
                        string_ptr_unordered_map_t* bcs, 
                        st_func_t* solution,
                        options_t* options)
{
  model_t* model = poisson_model_new(options);
  poisson_t* p = model_context(model);
  p->mesh = mesh;
  p->lambda = lambda;
  p->rhs = rhs;
  if (p->bcs != NULL)
    string_ptr_unordered_map_free(p->bcs);
  p->bcs = bcs;
  p->solution = solution;

  return model;
}

