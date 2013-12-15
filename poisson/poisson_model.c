// Copyright (c) 2012-2013, Jeffrey N. Johnson
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
#include "core/constant_st_func.h"
#include "core/boundary_cell_map.h"
#include "core/write_silo.h"
#include "geometry/interpreter_register_geometry_functions.h"
#include "integrators/nonlinear_solver.h"
#include "poisson/poisson_model.h"
#include "poisson/poisson_bc.h"
#include "poisson/interpreter_register_poisson_functions.h"
#include "poisson/register_poisson_benchmarks.h"

// Poisson model structure.
typedef struct 
{
  mesh_t* mesh;             
  point_cloud_t* point_cloud;
  adj_graph_t* graph;
  st_func_t* rhs;           // Right-hand side function.
  st_func_t* lambda;        // "Conduction" operator. 
  double* phi;              // Solution array.
  st_func_t* solution;      // Analytic solution (if non-NULL).

  string_ptr_unordered_map_t* bcs; // Boundary conditions.
  boundary_cell_map_t* boundary_cells; // Boundary cell info

  // Nonlinear solver that integrates Poisson's equation.
  nonlinear_solver_t* solver;
} poisson_t;

static void poisson_advance(void* context, double t, double dt)
{
  poisson_t* p = context;
  nonlinear_solver_solve(p->solver, t+dt, p->phi);
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
  p->lambda = interpreter_get_scalar_function(interp, "lambda");
  if (p->lambda == NULL)
  {
    double one = 1.0;
    p->lambda = constant_st_func_new(1, &one);
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

static int fv_poisson_residual(N_Vector u, N_Vector F, void* context)
{
  // FIXME
  return 0;
}

static int fvpm_poisson_residual(N_Vector u, N_Vector F, void* context)
{
  // FIXME
  return 0;
}

// Accesses the adjacency graph for the Poisson solver.
static adj_graph_t* get_graph(void* context)
{
  poisson_t* p = context;
  return p->graph;
}

static void poisson_init(void* context, double t)
{
  poisson_t* p = context;

  // If the model has been previously initialized, clean everything out.
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

  if (p->solver != NULL)
  {
    nonlinear_solver_free(p->solver);
    p->solver = NULL;
  }

  if (p->mesh != NULL)
  {
    // Initialize the solution vector.
    p->phi = malloc(sizeof(double)*p->mesh->num_cells);
    memset(p->phi, 0, sizeof(double)*p->mesh->num_cells);

    // Gather information about boundary cells.
    p->boundary_cells = boundary_cell_map_from_mesh_and_bcs(p->mesh, p->bcs);

    // Initialize the nonlinear solver.
    nonlinear_solver_vtable vtable = {.eval = fv_poisson_residual, .dtor = NULL, .graph = get_graph};
    p->solver = nonlinear_solver_new("Poisson (FV)", p, vtable, GMRES);
  }
  else
  {
    // Initialize the solution vector.
    p->phi = malloc(sizeof(double)*p->point_cloud->num_points);
    memset(p->phi, 0, sizeof(double)*p->point_cloud->num_points);

    // Gather information about boundary cells.
    //p->boundary_cells = boundary_cell_map_from_mesh_and_bcs(p->mesh, p->bcs);

    // Initialize the nonlinear solver.
    nonlinear_solver_vtable vtable = {.eval = fvpm_poisson_residual, .dtor = NULL, .graph = get_graph};
    p->solver = nonlinear_solver_new("Poisson (FVPM)", p, vtable, GMRES);
  }

  // Now we simply solve the problem for the initial time.
  nonlinear_solver_solve(p->solver, t, p->phi);

}

static void poisson_plot(void* context, const char* prefix, const char* directory, double t, int step)
{
  ASSERT(context != NULL);
  poisson_t* p = context;

  string_ptr_unordered_map_t* cell_fields = string_ptr_unordered_map_new();
  string_ptr_unordered_map_insert(cell_fields, "phi", p->phi);

  // If we are given an analytic solution, write it and the solution error.
  if (p->solution != NULL)
  {
    double *soln = malloc(sizeof(double) * p->mesh->num_cells), 
           *error = malloc(sizeof(double) * p->mesh->num_cells);
    for (int c = 0; c < p->mesh->num_cells; ++c)
    {
      st_func_eval(p->solution, &p->mesh->cell_centers[c], t, &soln[c]);
      error[c] = p->phi[c] - soln[c];
    }
    string_ptr_unordered_map_insert_with_v_dtor(cell_fields, "solution", soln, DTOR(free));
    string_ptr_unordered_map_insert_with_v_dtor(cell_fields, "error", error, DTOR(free));
  }
  if (p->mesh != NULL)
  {
    write_silo_mesh(p->mesh, NULL, NULL, NULL, cell_fields, 
                    prefix, directory, 0, 0.0, MPI_COMM_SELF, 1, 0);
  }
  else
  {
    write_silo_points(p->point_cloud->point_coords, 
                      p->point_cloud->num_points, 
                      cell_fields, prefix, directory, 0, 0.0, MPI_COMM_SELF, 1, 0);
  }
}

static void poisson_save(void* context, const char* prefix, const char* directory, double t, int step)
{
  ASSERT(context != NULL);
  poisson_t* p = context;

  string_ptr_unordered_map_t* cell_fields = string_ptr_unordered_map_new();
  string_ptr_unordered_map_insert(cell_fields, "phi", p->phi);
  if (p->mesh != NULL)
  {
    write_silo_mesh(p->mesh, NULL, NULL, NULL, cell_fields, 
                    prefix, directory, 0, 0.0, MPI_COMM_SELF, 1, 0);
  }
  else
  {
    write_silo_points(p->point_cloud->point_coords, 
                      p->point_cloud->num_points, 
                      cell_fields, prefix, directory, 0, 0.0, MPI_COMM_SELF, 1, 0);
  }
}

static void poisson_compute_error_norms(void* context, st_func_t* solution, double t, double* lp_norms)
{
  poisson_t* p = context;
  double Linf = 0.0, L1 = 0.0, L2 = 0.0;
  if (p->mesh != NULL)
  {
    for (int c = 0; c < p->mesh->num_cells; ++c)
    {
      double phi_sol;
      st_func_eval(solution, &p->mesh->cell_centers[c], t, &phi_sol);
      double V = p->mesh->cell_volumes[c];
      double err = fabs(p->phi[c] - phi_sol);
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
      double phi_sol;
      st_func_eval(solution, &p->point_cloud->point_coords[i], t, &phi_sol);
      double V = 1.0; // FIXME: What's the volume factor?
      double err = fabs(p->phi[i] - phi_sol);
      Linf = (Linf < err) ? err : Linf;
      L1 += err*V;
      L2 += err*err*V*V;
    }
  }

  L2 = sqrt(L2);
  lp_norms[0] = Linf;
  lp_norms[1] = L1;
  lp_norms[2] = L2;
}

static void poisson_dtor(void* context)
{
  poisson_t* p = context;

  // Destroy BC table.
  string_ptr_unordered_map_free(p->bcs);

  if (p->mesh != NULL)
    mesh_free(p->mesh);
  else if (p->point_cloud != NULL)
    point_cloud_free(p->point_cloud);

  if (p->phi != NULL)
    free(p->phi);

  if (p->boundary_cells != NULL)
    boundary_cell_map_free(p->boundary_cells);

  if (p->solver != NULL)
    nonlinear_solver_free(p->solver);

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
                                             {"lambda", INTERPRETER_SCALAR_FUNCTION},
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

