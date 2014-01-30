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
#include "core/point_cloud.h"
#include "core/least_squares.h"
#include "core/linear_algebra.h"
#include "core/write_silo.h"
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
  point_cloud_t* point_cloud;
  adj_graph_t* graph;
  st_func_t* rhs;           // Right-hand side function.
  st_func_t* lambda;        // "Conduction" operator (symmetric tensor). 
  real_t* phi;              // Solution array.
  st_func_t* solution;      // Analytic solution (if non-NULL).

  string_ptr_unordered_map_t* bcs; // Boundary conditions.
  boundary_cell_map_t* boundary_cells; // Boundary cell -> BC mapping

  // Here we store the fluxes between cells/points, and keep track of 
  // which fluxes have already been computed. The fluxes for each face of 
  // each cell are computed and stored in compressed row storage (CRS) format,
  // with each face having a distinct value for each of its cells.
  real_t* cell_face_fluxes;
  int* cell_face_flux_offsets;

  // Contributions from the source term (rhs).
  real_t* cell_sources;

  // Polynomial fit.
  polynomial_fit_t* poly_fit;

  // Quadrature rules -- regular and "special" (to accommodate symmetry).
  polyhedron_integrator_t* poly_quad_rule;
  int_ptr_unordered_map_t* special_quad_rules;

  // Nonlinear solver that integrates Poisson's equation.
  nonlinear_integrator_t* solver;

  // Number of iterations for most recent nonlinear solve.
  int num_iterations;
} poisson_t;

static void poisson_advance(void* context, real_t t, real_t dt)
{
  poisson_t* p = context;
  nonlinear_integrator_solve(p->solver, t+dt, p->phi, &p->num_iterations);
}

static void poisson_read_input(void* context, interpreter_t* interp, options_t* options)
{
  poisson_t* p = context;
  p->mesh = interpreter_get_mesh(interp, "mesh");
  if (p->mesh == NULL)
  {
    int N;
    point_t* points = interpreter_get_pointlist(interp, "points", &N);
    if (points == NULL)
      polymec_error("poisson: Neither points nor mesh were specified.");
    else
      p->point_cloud = point_cloud_new(MPI_COMM_WORLD, points, N);
  }
  p->rhs = interpreter_get_scalar_function(interp, "rhs");
  if (p->rhs == NULL)
    polymec_error("poisson: No right hand side (rhs) was specified.");
  p->lambda = interpreter_get_sym_tensor_function(interp, "lambda");
  if (p->lambda == NULL)
  {
    real_t ones[6] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    p->lambda = constant_st_func_new(6, ones);
  }
  p->bcs = interpreter_get_table(interp, "bcs");
  if (p->bcs == NULL)
    polymec_error("poisson: No table of boundary conditions (bcs) was specified.");

  if (p->mesh != NULL)
  {
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
}

static int fv_poisson_residual(void* context, real_t t, real_t* u, real_t* F)
{
  poisson_t* p = context;

  // Set up arrays for storing the nodes of faces attached to cells.
  static const int max_num_faces_per_cell = 30;
  static const int max_num_nodes_per_face = 30;
  point_t face_nodes[max_num_faces_per_cell * max_num_nodes_per_face];
  int face_node_offsets[max_num_faces_per_cell];

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
    {
      // Set the integration domain for this polyhedral cell.
      int fpos = 0, face, num_faces = 0, fn_offset = 0;
      while (mesh_cell_next_face(p->mesh, cell, &fpos, &face))
      {
        int npos = 0, node;
        face_node_offsets[face] = fn_offset;
        while (mesh_next_face_node(p->mesh, face, &npos, &node))
        {
          face_nodes[fn_offset] = p->mesh->nodes[node];
          ++fn_offset;
        }
        ++num_faces;
      }
      face_node_offsets[num_faces] = fn_offset;

      polyhedron_integrator_set_domain(quad_rule, face_nodes, 
                                       face_node_offsets, num_faces);
    }

    // Find a least-squares fit for the solution in the vicinity of this cell.
    polynomial_fit_compute(p->poly_fit, u, cell);

    // Face fluxes.
    int fpos = 0, face, offset = p->cell_face_flux_offsets[cell];
    while (mesh_cell_next_face(p->mesh, cell, &fpos, &face))
    {
      real_t face_flux = 0.0;

      // Go over the quadrature points in the face and compute the flux
      // at each point, accumulating the integral in face_flux.
      int qpos = 0;
      point_t xq;
      vector_t nq;
      real_t wq;
      while (polyhedron_integrator_next_surface_point(quad_rule, &qpos, &xq, &nq, &wq))
      {
        // Compute the gradient of the solution using the polynomial fit.
        vector_t grad_phi;
        polynomial_fit_eval_deriv(p->poly_fit, &xq, 1, 0, 0, &grad_phi.x);
        polynomial_fit_eval_deriv(p->poly_fit, &xq, 0, 1, 0, &grad_phi.y);
        polynomial_fit_eval_deriv(p->poly_fit, &xq, 0, 0, 1, &grad_phi.z);

        // Evaluate the conduction operator lambda.
        real_t lambda[6];
        st_func_eval(p->lambda, &xq, t, lambda);

        // Form the flux contribution.
        vector_t n_o_lambda = {.x = nq.x * lambda[0] + nq.y * lambda[3] + nq.z * lambda[2],
                               .y = nq.x * lambda[1] + nq.y * lambda[2] + nq.z * lambda[4],
                               .z = nq.x * lambda[2] + nq.y * lambda[4] + nq.z * lambda[5]};
        face_flux += wq * vector_dot(&n_o_lambda, &grad_phi);
      }

      p->cell_face_fluxes[offset] = face_flux;
      F[cell] -= face_flux;
      ++offset;
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
    int cpos = 0, other_cell, offset1 = p->cell_face_flux_offsets[cell];
    while (mesh_cell_next_neighbor(p->mesh, cell, &cpos, &other_cell))
    {
      if ((other_cell != -1) && (cell < other_cell))
      {
        // Retrieve the fluxes at face 1 and face 2.
        real_t flux1 = p->cell_face_fluxes[offset1];
        int f = offset1 - p->cell_face_flux_offsets[cell];
        int offset2 = p->cell_face_flux_offsets[other_cell] + f;
        real_t flux2 = p->cell_face_fluxes[offset2];

        // Either they are both zero, or they have opposite sign.
        ASSERT((flux1 == flux2 == 0.0) || (SIGN(flux1) == -SIGN(flux2)));

        // Now average their magnitudes and make sure that they agree so 
        // that our fluxes are conservative.
        real_t avg_flux = 0.5 * (fabs(flux1) + fabs(flux2));
        p->cell_face_fluxes[offset1] = SIGN(flux1) * avg_flux;
        p->cell_face_fluxes[offset2] = SIGN(flux2) * avg_flux;
        ++offset1;
      }
    }

    // Now compute the residual in this cell.
    F[cell] = p->cell_sources[cell];
    int fpos = 0, face, offset = p->cell_face_flux_offsets[cell];
    while (mesh_cell_next_face(p->mesh, cell, &fpos, &face))
    {
      F[cell] -= p->cell_face_fluxes[offset];
      ++offset;
    }
  }

  return 0;
}

static int fvpm_poisson_residual(void* context, real_t t, real_t* u, real_t* F)
{
  // FIXME
  return 0;
}

// This helper clears all the data structures that are tied to a discrete domain.
static void poisson_clear(poisson_t* p)
{
  if (p->cell_face_fluxes != NULL)
  {
    free(p->cell_face_fluxes);
    p->cell_face_fluxes = NULL;
  }
  if (p->cell_face_flux_offsets != NULL)
  {
    free(p->cell_face_flux_offsets);
    p->cell_face_flux_offsets = NULL;
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

  if (p->special_quad_rules != NULL)
  {
    int_ptr_unordered_map_free(p->special_quad_rules);
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

  if (p->mesh != NULL)
  {
    int N = p->mesh->num_cells;

    // Initialize the solution vector.
    p->phi = malloc(sizeof(real_t)*N);
    memset(p->phi, 0, sizeof(real_t)*N);

    // Gather information about boundary cells.
    p->boundary_cells = boundary_cell_map_from_mesh_and_bcs(p->mesh, p->bcs);

    // Initialize the nonlinear solver.
    nonlinear_integrator_vtable vtable = {.eval = fv_poisson_residual, 
                                          .dtor = NULL};
    p->solver = bicgstab_nonlinear_integrator_new("Poisson (FV)", p, MPI_COMM_WORLD, N, vtable, LINE_SEARCH, 15);

    // For now, Use LU preconditioning with the same residual function.
    // Use LU preconditioning with the same residual function.
    preconditioner_t* lu_precond = lu_preconditioner_new(p, fv_poisson_residual, p->graph);
    nonlinear_integrator_set_preconditioner(p->solver, lu_precond);

    // Allocate storage for cell face fluxes.
    p->cell_face_flux_offsets = malloc(sizeof(int)*(N+1));
    p->cell_face_flux_offsets[0] = 0;
    for (int c = 0; c < N; ++c)
      p->cell_face_flux_offsets[c+1] = p->cell_face_flux_offsets[c] + mesh_cell_num_faces(p->mesh, c);
    p->cell_face_fluxes = malloc(sizeof(real_t)*p->cell_face_flux_offsets[N]);
    memset(p->cell_face_fluxes, 0, sizeof(real_t)*p->cell_face_flux_offsets[N]);

    // Allocate storage for cell source contributions.
    p->cell_sources = malloc(sizeof(real_t)*N);
    memset(p->cell_sources, 0, sizeof(real_t)*N);
  }
  else
  {
    int N = p->point_cloud->num_points;

    // Initialize the solution vector.
    p->phi = malloc(sizeof(real_t)*N);
    memset(p->phi, 0, sizeof(real_t)*N);

    // Gather information about boundary cells.
    //p->boundary_cells = boundary_cell_map_from_mesh_and_bcs(p->mesh, p->bcs);

    // Initialize the nonlinear solver.
    nonlinear_integrator_vtable vtable = {.eval = fvpm_poisson_residual, 
                                          .dtor = NULL};
    p->solver = bicgstab_nonlinear_integrator_new("Poisson (FVPM)", p, MPI_COMM_WORLD, N, vtable, LINE_SEARCH, 15);

    // For now, Use LU preconditioning with the same residual function.
    preconditioner_t* lu_precond = lu_preconditioner_new(p, fvpm_poisson_residual, p->graph);
    nonlinear_integrator_set_preconditioner(p->solver, lu_precond);
  }

  // Now we simply solve the problem for the initial time.
  nonlinear_integrator_solve(p->solver, t, p->phi, &p->num_iterations);
}

static void poisson_plot(void* context, const char* prefix, const char* directory, real_t t, int step)
{
  ASSERT(context != NULL);
  poisson_t* p = context;

  string_ptr_unordered_map_t* cell_fields = string_ptr_unordered_map_new();
  string_ptr_unordered_map_insert(cell_fields, "phi", p->phi);

  // If we are given an analytic solution, write it and the solution error.
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
  if (p->mesh != NULL)
    write_silo_mesh(p->mesh, cell_fields, prefix, directory, 0, 0.0, MPI_COMM_SELF, 1, 0);
  else
  {
    write_silo_points(p->point_cloud->point_coords, 
                      p->point_cloud->num_points, 
                      cell_fields, prefix, directory, 0, 0.0, MPI_COMM_SELF, 1, 0);
  }
}

static void poisson_save(void* context, const char* prefix, const char* directory, real_t t, int step)
{
  ASSERT(context != NULL);
  poisson_t* p = context;

  string_ptr_unordered_map_t* cell_fields = string_ptr_unordered_map_new();
  string_ptr_unordered_map_insert(cell_fields, "phi", p->phi);
  if (p->mesh != NULL)
    write_silo_mesh(p->mesh, cell_fields, prefix, directory, 0, 0.0, MPI_COMM_SELF, 1, 0);
  else
  {
    write_silo_points(p->point_cloud->point_coords, 
                      p->point_cloud->num_points, 
                      cell_fields, prefix, directory, 0, 0.0, MPI_COMM_SELF, 1, 0);
  }
}

static void poisson_compute_error_norms(void* context, st_func_t* solution, real_t t, real_t* lp_norms)
{
  poisson_t* p = context;
  real_t Linf = 0.0, L1 = 0.0, L2 = 0.0;
  if (p->mesh != NULL)
  {
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
  }
  else
  {
    for (int i = 0; i < p->point_cloud->num_points; ++i)
    {
      real_t phi_sol;
      st_func_eval(solution, &p->point_cloud->point_coords[i], t, &phi_sol);
      real_t V = 1.0; // FIXME: What's the volume factor?
      real_t err = fabs(p->phi[i] - phi_sol);
      Linf = (Linf < err) ? err : Linf;
      L1 += err*V;
      L2 += err*err*V*V;
    }
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
  else if (p->point_cloud != NULL)
    point_cloud_free(p->point_cloud);
  if (p->poly_fit != NULL)
    polynomial_fit_free(p->poly_fit);

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
  p->rhs = NULL;
  p->phi = NULL;
  p->solver = NULL;
  p->bcs = string_ptr_unordered_map_new();
  p->solution = NULL;

  p->boundary_cells = boundary_cell_map_new();
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

model_t* create_fv_poisson(mesh_t* mesh,
                           st_func_t* lambda,
                           st_func_t* rhs,
                           string_ptr_unordered_map_t* bcs, 
                           st_func_t* solution,
                           options_t* options)
{
  model_t* model = poisson_model_new(options);
  poisson_t* pm = model_context(model);
  pm->mesh = mesh;
  pm->lambda = lambda;
  pm->rhs = rhs;
  if (pm->bcs != NULL)
    string_ptr_unordered_map_free(pm->bcs);
  pm->bcs = bcs;
  pm->solution = solution;

  return model;
}

model_t* create_fvpm_poisson(point_cloud_t* point_cloud,
                             st_func_t* lambda,
                             st_func_t* rhs,
                             string_ptr_unordered_map_t* bcs, 
                             st_func_t* solution,
                             options_t* options)
{
  model_t* model = poisson_model_new(options);
  poisson_t* pm = model_context(model);
  pm->point_cloud = point_cloud;
  pm->lambda = lambda;
  pm->rhs = rhs;
  if (pm->bcs != NULL)
    string_ptr_unordered_map_free(pm->bcs);
  pm->bcs = bcs;
  pm->solution = solution;

  return model;
}

