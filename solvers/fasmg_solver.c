// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/linear_algebra.h"
#include "solvers/fasmg_solver.h"

struct fasmg_solver_t 
{
  fasmg_operator_t* A;
  fasmg_coarsener_t* coarsener;
  fasmg_restrictor_t* restrictor;
  fasmg_prolongator_t* prolongator;
  fasmg_cycle_t* cycle;

  int max_cycles;
  real_t max_res_norm;
};

struct fasmg_grid_t 
{
  void* data;
  size_t num_dof;
  fasmg_grid_vtable vtable;
  bool free_data;

  fasmg_grid_t* coarser;
  fasmg_grid_t* finer;
};

struct fasmg_operator_t 
{
  char* name;
  void* context;
  fasmg_operator_vtable vtable;
};

struct fasmg_coarsener_t 
{
  char* name;
  void* context;
  fasmg_coarsener_vtable vtable;
};

struct fasmg_restrictor_t 
{
  char* name;
  void* context;
  fasmg_restrictor_vtable vtable;
};

struct fasmg_prolongator_t 
{
  char* name;
  void* context;
  fasmg_prolongator_vtable vtable;
};

struct fasmg_cycle_t
{
  char* name;
  void* context;
  fasmg_cycle_vtable vtable;
};

fasmg_solver_t* fasmg_solver_new(fasmg_operator_t* A,
                                 fasmg_coarsener_t* coarsener,
                                 fasmg_restrictor_t* restrictor,
                                 fasmg_prolongator_t* prolongator,
                                 fasmg_cycle_t* cycle)
{
  ASSERT(A != NULL);
  ASSERT(coarsener != NULL);
  ASSERT(restrictor != NULL);
  ASSERT(prolongator != NULL);
  ASSERT(cycle != NULL);

  fasmg_solver_t* solver = polymec_malloc(sizeof(fasmg_solver_t));
  solver->A = A;
  solver->coarsener = coarsener;
  solver->restrictor = restrictor;
  solver->prolongator = prolongator;
  solver->cycle = cycle;
  solver->max_cycles = 30;
  solver->max_res_norm = 1e-10;
  return solver;
}

void fasmg_solver_free(fasmg_solver_t* solver)
{
  fasmg_cycle_free(solver->cycle);
  fasmg_prolongator_free(solver->prolongator);
  fasmg_restrictor_free(solver->restrictor);
  fasmg_coarsener_free(solver->coarsener);
  fasmg_operator_free(solver->A);
  polymec_free(solver);
}

void fasmg_solver_set_max_cycles(fasmg_solver_t* solver,
                                 int max_cycles)
{
  ASSERT(max_cycles >= 1);
  solver->max_cycles = max_cycles;
}

void fasmg_solver_set_max_residual_norm(fasmg_solver_t* solver,
                                        real_t max_residual_norm)
{
  ASSERT(max_residual_norm > 0.0);
  solver->max_res_norm = max_residual_norm;
}

static fasmg_grid_t* fasmg_grid_new(void* data,
                                    size_t num_dof,
                                    fasmg_grid_vtable vtable)
{
  ASSERT(num_dof > 0);

  fasmg_grid_t* grid = polymec_malloc(sizeof(fasmg_grid_t));
  grid->data = data;
  grid->num_dof = num_dof;
  grid->vtable = vtable;
  grid->free_data = true;
  grid->coarser = grid->finer = NULL;
  return grid;
}

fasmg_grid_t* fasmg_solver_grid(fasmg_solver_t* solver,
                                void* discretization,
                                size_t num_dof,
                                fasmg_grid_vtable vtable)
{
  ASSERT(vtable.l2_norm != NULL);
  fasmg_grid_t* grid = fasmg_grid_new(discretization, num_dof, vtable); 
  grid->free_data = false;
  return grid;
}

bool fasmg_solver_solve(fasmg_solver_t* solver,
                        fasmg_grid_t* grid,
                        real_t* B,
                        real_t* X,
                        real_t* residual_norm,
                        int* num_cycles)
{
  int n_cycles = 0;

  // Compute the initial residual norm.
  size_t N = grid->num_dof;
  real_t R[N];
  fasmg_operator_compute_residual(solver->A, grid, B, X, R);
  *residual_norm = fasmg_grid_l2_norm(grid, R);

  while ((n_cycles < solver->max_cycles) && (*residual_norm > solver->max_res_norm))
  {
    real_t prev_res_norm = *residual_norm;

    // Do a cycle (FMG first, then "regular" cycles).
    if (n_cycles == 0)
      fasmg_solver_fmg(solver, 1, grid, B, X);
    else
      fasmg_solver_cycle(solver, grid, B, X);
    ++n_cycles;

    // Compute the residual norm.
    fasmg_operator_compute_residual(solver->A, grid, B, X, R);
    *residual_norm = fasmg_grid_l2_norm(grid, R);
    real_t conv_ratio = *residual_norm / prev_res_norm;

    log_detail("fasmg_solver_solve: cycle %d, residual norm = %g, conv. ratio = %g", 
               n_cycles, *residual_norm, conv_ratio);
  }
  *num_cycles = n_cycles;
  return ((n_cycles <= solver->max_cycles) && 
          (*residual_norm < solver->max_res_norm));
}

void fasmg_solver_cycle(fasmg_solver_t* solver,
                        fasmg_grid_t* grid,
                        real_t* B,
                        real_t* X)
{
  // Create a coarse hierarchy of grids for this one if needed.
  fasmg_coarsener_coarsen(solver->coarsener, grid);

  // Cycle.
  fasmg_cycle_execute(solver->cycle, solver->A, grid, 
                      solver->prolongator, solver->restrictor, B, X);
}

static void execute_fmg(fasmg_operator_t* A, 
                        fasmg_grid_t* grid, 
                        fasmg_prolongator_t* prolongator,
                        fasmg_restrictor_t* restrictor,
                        int nu_0,
                        fasmg_cycle_t* cycle, 
                        real_t* B,
                        real_t* X)
{
  size_t N = grid->num_dof;

  // This follows the recipe in p. 43 of the 2nd edition of A Multigrid 
  // Tutorial by Briggs, Henson, and McCormick (SIAM 1999).
  if (grid->coarser == NULL)
    memset(X, 0, sizeof(real_t) * N);
  else
  {
    // Restrict B to our coarser grid.
    size_t N_coarse = grid->coarser->num_dof;
    real_t B_coarse[N_coarse];
    fasmg_restrictor_project(restrictor, grid, B, B_coarse);

    // Cycle on the coarse solution.
    real_t X_coarse[N_coarse];
    log_indent(LOG_DEBUG);
    execute_fmg(A, grid->coarser, prolongator, restrictor, nu_0, cycle,
                B_coarse, X_coarse);
    log_unindent(LOG_DEBUG);

    // Now prolongate our coarse solution to the present one.
    fasmg_prolongator_interpolate(prolongator, grid->coarser, X_coarse, X);
  }

  // V cycle nu_0 times.
  for (int i = 0; i < nu_0; ++i)
    fasmg_cycle_execute(cycle, A, grid, prolongator, restrictor, B, X);
}

void fasmg_solver_fmg(fasmg_solver_t* solver,
                      int nu_0,
                      fasmg_grid_t* grid,
                      real_t* B,
                      real_t* X)
{
  // Create a coarse hierarchy of grids for this one if needed.
  fasmg_coarsener_coarsen(solver->coarsener, grid);

  // Do the full multigrid step.
  execute_fmg(solver->A, grid, solver->prolongator, solver->restrictor, 
              nu_0, solver->cycle, B, X);
}

fasmg_operator_t* fasmg_solver_operator(fasmg_solver_t* solver)
{
  return solver->A;
}

fasmg_coarsener_t* fasmg_solver_coarsener(fasmg_solver_t* solver)
{
  return solver->coarsener;
}

fasmg_restrictor_t* fasmg_solver_restrictor(fasmg_solver_t* solver)
{
  return solver->restrictor;
}

fasmg_prolongator_t* fasmg_solver_prolongator(fasmg_solver_t* solver)
{
  return solver->prolongator;
}

fasmg_cycle_t* fasmg_solver_cycler(fasmg_solver_t* solver)
{
  return solver->cycle;
}

void fasmg_grid_free(fasmg_grid_t* grid)
{
  // Recursively delete the coarser grids.
  if (grid->coarser != NULL)
    fasmg_grid_free(grid->coarser);

  // Delete this one.
  if ((grid->free_data) && 
      (grid->vtable.dtor != NULL) && 
      (grid->data != NULL))
    grid->vtable.dtor(grid->data);
  polymec_free(grid);
}

size_t fasmg_grid_num_dof(fasmg_grid_t* grid)
{
  return grid->num_dof;
}

real_t fasmg_grid_l2_norm(fasmg_grid_t* grid, real_t* V)
{
  return grid->vtable.l2_norm(grid->data, V);
}

fasmg_grid_t* fasmg_grid_coarser(fasmg_grid_t* grid)
{
  return grid->coarser;
}

fasmg_grid_t* fasmg_grid_finer(fasmg_grid_t* grid)
{
  return grid->finer;
}

fasmg_operator_t* fasmg_operator_new(const char* name, 
                                     void* context,
                                     fasmg_operator_vtable vtable)
{
  fasmg_operator_t* A = polymec_malloc(sizeof(fasmg_operator_t));
  A->name = string_dup(name);
  A->context = context;
  A->vtable = vtable;
  return A;
}

void fasmg_operator_free(fasmg_operator_t* A)
{
  if ((A->vtable.dtor != NULL) && (A->context != NULL))
    A->vtable.dtor(A->context);
  string_free(A->name);
  polymec_free(A);
}

char* fasmg_operator_name(fasmg_operator_t* A)
{
  return A->name;
}

void* fasmg_operator_context(fasmg_operator_t* A)
{
  return A->context;
}

void fasmg_operator_apply(fasmg_operator_t* A,
                          fasmg_grid_t* grid,
                          real_t* X, 
                          real_t* AX)
{
  A->vtable.apply(A->context, grid->data, X, AX);
}

void fasmg_operator_relax(fasmg_operator_t* A,
                          fasmg_grid_t* grid,
                          real_t* B,
                          real_t* X)
{
  A->vtable.relax(A->context, grid->data, B, X);
}

bool fasmg_operator_has_direct_solve(fasmg_operator_t* A)
{
  return (A->vtable.solve_directly != NULL);
}

void fasmg_operator_solve_directly(fasmg_operator_t* A,
                                   fasmg_grid_t* grid,
                                   real_t* B,
                                   real_t* X)
{
  ASSERT(A->vtable.solve_directly != NULL);
  A->vtable.solve_directly(A->context, grid->data, B, X);
}

void fasmg_operator_compute_residual(fasmg_operator_t* A,
                                     fasmg_grid_t* grid,
                                     real_t* B,
                                     real_t* X,
                                     real_t* R)
{
  A->vtable.apply(A->context, grid->data, X, R); // R <- AX
  for (size_t i = 0; i < grid->num_dof; ++i)
    R[i] = B[i] - R[i]; // R <- B - AX
}

fasmg_coarsener_t* fasmg_coarsener_new(const char* name, 
                                       void* context,
                                       fasmg_coarsener_vtable vtable)
{
  fasmg_coarsener_t* coarsener = polymec_malloc(sizeof(fasmg_coarsener_t));
  coarsener->name = string_dup(name);
  coarsener->context = context;
  coarsener->vtable = vtable;
  return coarsener;
}

void fasmg_coarsener_free(fasmg_coarsener_t* coarsener)
{
  if ((coarsener->vtable.dtor != NULL) && (coarsener->context != NULL))
    coarsener->vtable.dtor(coarsener->context);
  string_free(coarsener->name);
  polymec_free(coarsener);
}

char* fasmg_coarsener_name(fasmg_coarsener_t* coarsener)
{
  return coarsener->name;
}

void* fasmg_coarsener_context(fasmg_coarsener_t* coarsener)
{
  return coarsener->context;
}

bool fasmg_coarsener_can_coarsen(fasmg_coarsener_t* coarsener,
                                 fasmg_grid_t* grid)
{
  return coarsener->vtable.can_coarsen(coarsener->context, grid->data);
}

void fasmg_coarsener_coarsen(fasmg_coarsener_t* coarsener,
                             fasmg_grid_t* grid)
{
  // Create a hierarachy of coarse grids if we haven't already.
  fasmg_grid_t* G = grid;
  while ((G->coarser == NULL) && 
         fasmg_coarsener_can_coarsen(coarsener, G))
  {
    size_t coarse_dof = 0;
    void* coarse_data = coarsener->vtable.coarser_grid(coarsener->context, G->data, &coarse_dof);
    fasmg_grid_t* coarse_grid = fasmg_grid_new(coarse_data, coarse_dof, G->vtable);
    G->coarser = coarse_grid;
    coarse_grid->finer = G;
    G = coarse_grid;
  }
}

fasmg_restrictor_t* fasmg_restrictor_new(const char* name, 
                                         void* context,
                                         fasmg_restrictor_vtable vtable)
{
  fasmg_restrictor_t* restrictor = polymec_malloc(sizeof(fasmg_restrictor_t));
  restrictor->name = string_dup(name);
  restrictor->context = context;
  restrictor->vtable = vtable;
  return restrictor;
}

void fasmg_restrictor_free(fasmg_restrictor_t* restrictor)
{
  if ((restrictor->vtable.dtor != NULL) && (restrictor->context != NULL))
    restrictor->vtable.dtor(restrictor->context);
  string_free(restrictor->name);
  polymec_free(restrictor);
}

char* fasmg_restrictor_name(fasmg_restrictor_t* restrictor)
{
  return restrictor->name;
}

void* fasmg_restrictor_context(fasmg_restrictor_t* restrictor)
{
  return restrictor->context;
}

void fasmg_restrictor_project(fasmg_restrictor_t* restrictor,
                              fasmg_grid_t* fine_grid,
                              real_t* fine_X,
                              real_t* coarse_X)
{
  restrictor->vtable.project(restrictor->context, fine_grid->data, 
                             fine_X, fine_grid->coarser->data, coarse_X);
}

fasmg_prolongator_t* fasmg_prolongator_new(const char* name, 
                                           void* context,
                                           fasmg_prolongator_vtable vtable)
{
  fasmg_prolongator_t* prolongator = polymec_malloc(sizeof(fasmg_prolongator_t));
  prolongator->name = string_dup(name);
  prolongator->context = context;
  prolongator->vtable = vtable;
  return prolongator;
}

void fasmg_prolongator_free(fasmg_prolongator_t* prolongator)
{
  if ((prolongator->vtable.dtor != NULL) && (prolongator->context != NULL))
    prolongator->vtable.dtor(prolongator->context);
  string_free(prolongator->name);
  polymec_free(prolongator);
}

char* fasmg_prolongator_name(fasmg_prolongator_t* prolongator)
{
  return prolongator->name;
}

void* fasmg_prolongator_context(fasmg_prolongator_t* prolongator)
{
  return prolongator->context;
}

void fasmg_prolongator_interpolate(fasmg_prolongator_t* prolongator,
                                   fasmg_grid_t* coarse_grid,
                                   real_t* coarse_X,
                                   real_t* fine_X)
{
  prolongator->vtable.interpolate(prolongator->context, coarse_grid->data, 
                                  coarse_X, coarse_grid->finer->data, fine_X);
}

fasmg_cycle_t* fasmg_cycle_new(const char* name, 
                               void* context,
                               fasmg_cycle_vtable vtable)
{
  fasmg_cycle_t* cycle = polymec_malloc(sizeof(fasmg_cycle_t));
  cycle->name = string_dup(name);
  cycle->context = context;
  cycle->vtable = vtable;
  return cycle;
}

void fasmg_cycle_free(fasmg_cycle_t* cycle)
{
  if ((cycle->vtable.dtor != NULL) && (cycle->context != NULL))
    cycle->vtable.dtor(cycle->context);
  string_free(cycle->name);
  polymec_free(cycle);
}

char* fasmg_cycle_name(fasmg_cycle_t* cycle)
{
  return cycle->name;
}

void* fasmg_cycle_context(fasmg_cycle_t* cycle)
{
  return cycle->context;
}

void fasmg_cycle_execute(fasmg_cycle_t* cycle,
                         fasmg_operator_t* A,
                         fasmg_grid_t* grid,
                         fasmg_prolongator_t* prolongator,
                         fasmg_restrictor_t* restrictor,
                         real_t* B,
                         real_t* X)
{
  cycle->vtable.execute(cycle->context, A, grid, prolongator, restrictor, B, X);
}

fasmg_cycle_t* v_fasmg_cycle_new(int nu_1, int nu_2)
{
  return mu_fasmg_cycle_new(1, nu_1, nu_2);
}

fasmg_cycle_t* w_fasmg_cycle_new(int nu_1, int nu_2)
{
  return mu_fasmg_cycle_new(2, nu_1, nu_2);
}

typedef struct
{
  int mu, nu_1, nu_2;
} mu_cycle_t;

static void mu_execute(void* context, 
                       fasmg_operator_t* A, 
                       fasmg_grid_t* grid, 
                       fasmg_prolongator_t* prolongator,
                       fasmg_restrictor_t* restrictor,
                       real_t* B,
                       real_t* X)
{
  mu_cycle_t* mu_cycle = context;

  // Relax or directly solve.
  if ((grid->coarser == NULL) && fasmg_operator_has_direct_solve(A))
  {
    log_debug("mu_cycle_execute: solving directly (N=%zd).", fasmg_grid_num_dof(grid));
    fasmg_operator_solve_directly(A, grid, B, X);
  }
  else
  {
    // Relax X nu_1 times.
    log_debug("mu_cycle_execute: relaxing %d times.", mu_cycle->nu_1);
    for (int i = 0; i < mu_cycle->nu_1; ++i)
      fasmg_operator_relax(A, grid, B, X);
  }

  if (grid->coarser != NULL)
  {
    // Compute the residual on this grid.
    size_t N = grid->num_dof;
    real_t R[N];
    fasmg_operator_compute_residual(A, grid, B, X, R);
    if (log_level() == LOG_DEBUG)
    {
      log_debug("mu_cycle_execute: residual norm is %g (N=%d)", 
                fasmg_grid_l2_norm(grid, R), N);
    }

    // Restrict the residual and our current solution to our coarser grid.
    size_t N_coarse = grid->coarser->num_dof;
    real_t R_coarse[N_coarse], X_coarse[N_coarse];
    fasmg_restrictor_project(restrictor, grid, R, R_coarse);
    fasmg_restrictor_project(restrictor, grid, X, X_coarse);

    // Form our right hand side on the coarse grid.
    real_t AX_coarse[N_coarse];
    fasmg_operator_apply(A, grid->coarser, X_coarse, AX_coarse);
    real_t B_coarse[N_coarse];
    for (size_t i = 0; i < N_coarse; ++i)
      B_coarse[i] = AX_coarse[i] + R_coarse[i];

    // Cycle on the coarse solution mu times to produce a solution Y_coarse.
    real_t Y_coarse[N_coarse];
    memcpy(Y_coarse, X_coarse, sizeof(real_t) * N_coarse);
    log_indent(LOG_DEBUG);
    for (int i = 0; i < mu_cycle->mu; ++i)
    {
      mu_execute(context, A, grid->coarser, prolongator, restrictor,
                 B_coarse, Y_coarse);
    }
    log_unindent(LOG_DEBUG);

    // Compute the coarse error approximation E_coarse = Y_coarse - X_coarse.
    real_t E_coarse[N_coarse];
    for (size_t i = 0; i < N_coarse; ++i)
      E_coarse[i] = Y_coarse[i] - X_coarse[i];

    // Prolongate the coarse error to our original grid.
    real_t E[N];
    fasmg_prolongator_interpolate(prolongator, grid->coarser, E_coarse, E);

    // Correct our current approximation.
    for (size_t i = 0; i < N; ++i)
      X[i] += E[i];
  }

  // Relax X nu_2 times.
  if ((grid->coarser != NULL) || !fasmg_operator_has_direct_solve(A))
  {
    log_debug("mu_cycle_execute: relaxing %d times.", mu_cycle->nu_2);
    for (int i = 0; i < mu_cycle->nu_2; ++i)
      fasmg_operator_relax(A, grid, B, X);
  }
}

fasmg_cycle_t* mu_fasmg_cycle_new(int mu, int nu_1, int nu_2)
{
  ASSERT(mu >= 0);
  ASSERT(nu_1 >= 0);
  ASSERT(nu_2 >= 0);

  mu_cycle_t* mu_cycle = polymec_malloc(sizeof(mu_cycle_t));
  mu_cycle->mu = mu;
  mu_cycle->nu_1 = nu_1;
  mu_cycle->nu_2 = nu_2;
  fasmg_cycle_vtable vtable = {.execute = mu_execute, .dtor = polymec_free};
  char name[1025];
  if (mu == 1)
    snprintf(name, 1024, "V cycle (nu_1 = %d, nu_2 = %d)", nu_1, nu_2);
  else if (mu == 2)
    snprintf(name, 1024, "W cycle (nu_1 = %d, nu_2 = %d)", nu_1, nu_2);
  else
  {
    snprintf(name, 1024, "mu cycle (mu = %d, nu_1 = %d, nu_2 = %d)", 
             mu, nu_1, nu_2);
  }
  return fasmg_cycle_new(name, mu_cycle, vtable);
}

