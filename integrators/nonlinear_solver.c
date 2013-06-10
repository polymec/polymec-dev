#include <float.h>
#include "integrators/nonlinear_solver.h"
#include "core/hypre_helpers.h"

struct nonlinear_function_t 
{
  void* context;
  char* name;
  nonlinear_function_vtable vtable;
};

struct nonlinear_timestepper_t
{
  void* context;
  char* name;
  int max_history_length, history_length;
  double* dt_history;
  double* error_history;
  int* iteration_history;
  double last_unsuccessful_dt;
  int num_failures;
  char* explanation;
  nonlinear_timestepper_vtable vtable;
};

struct nonlinear_solver_t 
{
  nonlinear_function_t* F;
  nonlinear_timestepper_t* timestepper;

  // Adjacency graph and coloring.
  adj_graph_t* graph;
  adj_graph_coloring_t* coloring;
  // Mapping of vertices to colors.
  int* colors;
  // Mapping of (color, row) to column.
  int** columns;

  int num_sites;

  // Linear system and solver.
  HYPRE_IJMatrix A;
  HYPRE_IJVector x, b;
  HYPRE_Solver solver;
  HYPRE_Solver precond;

  // Flag is set to true if the linear system above is initialized.
  bool initialized;

  // Index space.
  index_space_t* index_space;

  // Jacobian and residual coefficients.
  double_table_t* Jij;
  double *R, *deltaX;

  // Vectors for computing the Jacobian components on colored vertices. 
  double **Xps, **Rps, **ds, **Jds; 

  // Increment for approximating the Jacobian with finite differences.
  double delta;

  // Residual norm tolerance.
  double tolerance;

  // Maximum number of iterations before abject failure.
  int max_iters;
};

nonlinear_function_t* nonlinear_function_new(const char* name, 
                                             void* context,
                                             nonlinear_function_vtable vtable)
{
  ASSERT(vtable.eval_residual != NULL);

  nonlinear_function_t* F = malloc(sizeof(nonlinear_function_t));
  F->name = strdup(name);
  F->context = context;
  F->vtable = vtable;

  return F;
}

void nonlinear_function_free(nonlinear_function_t* F)
{
  if ((F->context != NULL) && (F->vtable.dtor != NULL))
    F->vtable.dtor(F->context);
  free(F->name);
  free(F);
}

char* nonlinear_function_name(nonlinear_function_t* F)
{
  return F->name;
}

void* nonlinear_function_context(nonlinear_function_t* F)
{
  return F->context;
}

void nonlinear_function_eval(nonlinear_function_t* F,
                             double t,
                             double* X,
                             double* R)
{
  ASSERT(X != NULL);
  ASSERT(R != NULL);
  F->vtable.eval_residual(F->context, t, X, R);
}

nonlinear_timestepper_t* nonlinear_timestepper_new(const char* name, 
                                                   void* context,
                                                   int max_history_length,
                                                   nonlinear_timestepper_vtable vtable)
{
  ASSERT(max_history_length >= 1);
  ASSERT(vtable.compute_dt != NULL);
  ASSERT(vtable.recompute_J != NULL);

  nonlinear_timestepper_t* timestepper = malloc(sizeof(nonlinear_timestepper_t));
  timestepper->name = strdup(name);
  timestepper->context = context;
  timestepper->history_length = 0;
  timestepper->max_history_length = max_history_length;
  timestepper->dt_history = malloc(sizeof(double) * max_history_length);
  timestepper->error_history = malloc(sizeof(double) * max_history_length);
  timestepper->iteration_history = malloc(sizeof(int) * max_history_length);
  timestepper->last_unsuccessful_dt = FLT_MAX;
  timestepper->num_failures = 0;
  timestepper->vtable = vtable;
  timestepper->explanation = malloc(sizeof(char) * 2048);
  
  return timestepper;
}

void nonlinear_timestepper_free(nonlinear_timestepper_t* timestepper)
{
  if ((timestepper->context != NULL) && (timestepper->vtable.dtor != NULL))
    timestepper->vtable.dtor(timestepper->context);
  free(timestepper->dt_history);
  free(timestepper->error_history);
  free(timestepper->iteration_history);
  free(timestepper->explanation);
  free(timestepper->name);
}

char* nonlinear_timestepper_name(nonlinear_timestepper_t* timestepper)
{
  return timestepper->name;
}

int nonlinear_timestepper_max_history_length(nonlinear_timestepper_t* timestepper)
{
  return timestepper->max_history_length;
}

void nonlinear_timestepper_clear_history(nonlinear_timestepper_t* timestepper)
{
  timestepper->history_length = 0;
}

double nonlinear_timestepper_step_size(nonlinear_timestepper_t* timestepper,
                                       char** explanation)
{
  ASSERT(explanation != NULL);
  snprintf(timestepper->explanation, 2048, "no explanation given");
  double dt = timestepper->vtable.compute_dt(timestepper->context, 
                                             timestepper->last_unsuccessful_dt,
                                             timestepper->num_failures,
                                             timestepper->dt_history + timestepper->history_length,
                                             timestepper->error_history + timestepper->history_length,
                                             timestepper->iteration_history + timestepper->history_length,
                                             timestepper->history_length,
                                             timestepper->explanation);
  *explanation = timestepper->explanation;
  return dt;
}

bool nonlinear_timestepper_recompute_jacobian(nonlinear_timestepper_t* timestepper)
{
  return timestepper->vtable.recompute_J(timestepper->context, 
                                         timestepper->dt_history + timestepper->history_length,
                                         timestepper->error_history + timestepper->history_length,
                                         timestepper->iteration_history + timestepper->history_length,
                                         timestepper->history_length);
}

void nonlinear_timestepper_converged(nonlinear_timestepper_t* timestepper,
                                     double step_size,
                                     double error_norm,
                                     int num_iterations)
{
  ASSERT(step_size > 0.0);
  ASSERT(error_norm >= 0.0);
  ASSERT(num_iterations > 0);

  if (timestepper->history_length == (timestepper->max_history_length - 1))
  {
    // We shift the history downward.
    for (int i = 1; i < timestepper->max_history_length; ++i)
    {
      timestepper->dt_history[i-1] = timestepper->dt_history[i];
      timestepper->error_history[i-1] = timestepper->error_history[i];
      timestepper->iteration_history[i-1] = timestepper->iteration_history[i];
    }
  }

  // Append the latest information to the history.
  timestepper->dt_history[timestepper->history_length] = step_size;
  timestepper->error_history[timestepper->history_length] = error_norm;
  timestepper->iteration_history[timestepper->history_length] = num_iterations;

  if (timestepper->history_length < (timestepper->max_history_length - 1))
    ++timestepper->history_length;

  // Reset the failure counter.
  timestepper->num_failures = 0;
}

void nonlinear_timestepper_failed(nonlinear_timestepper_t* timestepper,
                                  double step_size)
{
  // Record the time step size and increment the failure counter.
  timestepper->last_unsuccessful_dt = step_size;
  ++timestepper->num_failures;
}

nonlinear_solver_t* nonlinear_solver_new(nonlinear_function_t* F,
                                         nonlinear_timestepper_t* timestepper,
                                         adj_graph_t* graph)
{
  ASSERT(F != NULL);
  ASSERT(graph != NULL);

  nonlinear_solver_t* solver = malloc(sizeof(nonlinear_solver_t));
  solver->F = F;
  solver->timestepper = timestepper;
  solver->graph = graph;
  solver->coloring = NULL;
  solver->colors = NULL;
  solver->ds = NULL;
  solver->Jds = NULL;
  solver->Xps = NULL;
  solver->Rps = NULL;
  solver->columns = NULL;
  solver->num_sites = adj_graph_num_vertices(graph);
  solver->R = malloc(sizeof(double) * solver->num_sites);
  solver->deltaX = malloc(sizeof(double) * solver->num_sites);
  solver->index_space = index_space_from_low_and_high(adj_graph_comm(graph),
                                                      adj_graph_first_vertex(graph),
                                                      adj_graph_last_vertex(graph)+1);
  solver->initialized = false;
  solver->delta = FLT_EPSILON;
  solver->tolerance = 1e-8;
  solver->max_iters = 100;

  // Now set up the graph coloring stuff for the calculation of Jacobians.
  solver->coloring = adj_graph_coloring_new(solver->graph, SMALLEST_LAST); 
  int num_colors = adj_graph_coloring_num_colors(solver->coloring);
  log_detail("nonlinear_solver: graph coloring produced %d colors.", 
             num_colors);

  // Now that we have the coloring, we can set up our vectors {d} that 
  // allow us to compute the different colored portions of the Jacobian
  // in parallel.
  int num_vertices = solver->num_sites;
  solver->colors = malloc(num_vertices * sizeof(int));
  solver->columns = malloc(num_vertices * sizeof(int*));
  solver->Xps = malloc(num_colors * sizeof(double*));
  solver->Rps = malloc(num_colors * sizeof(double*));
  solver->ds = malloc(num_colors * sizeof(double*));
  solver->Jds = malloc(num_colors * sizeof(double*));
  for (int color = 0; color < num_colors; ++color)
  {
    solver->columns[color] = malloc(sizeof(int) * num_vertices);
    for (int i = 0; i < num_vertices; ++i)
      solver->columns[color][i] = -1;
    solver->Xps[color] = malloc(sizeof(double) * num_vertices);
    solver->Rps[color] = malloc(sizeof(double) * num_vertices);
    solver->ds[color] = malloc(sizeof(double) * num_vertices);
    solver->Jds[color] = malloc(sizeof(double) * num_vertices);
    memset(solver->Xps[color], 0, sizeof(double) * num_vertices);
    memset(solver->Rps[color], 0, sizeof(double) * num_vertices);
    memset(solver->ds[color], 0, sizeof(double) * num_vertices);
    memset(solver->Jds[color], 0, sizeof(double) * num_vertices);
    int pos = 0, vertex;
    while (adj_graph_coloring_next_vertex(solver->coloring, color, &pos, &vertex))
    {
      solver->colors[vertex] = color;
      solver->ds[color][vertex] = 1.0;
    }
  }

  // Map (color, row) pairs to columns using the graph.
  for (int i = 0; i < num_vertices; ++i)
  {
    int icolor = solver->colors[i];
    solver->columns[icolor][i] = i;
    int num_edges = adj_graph_num_edges(solver->graph, i);
    int *edges = adj_graph_edges(solver->graph, i);
    for (int e = 0; e < num_edges; ++e)
    {
      int j = edges[e];
      int jcolor = solver->colors[j];
      solver->columns[jcolor][i] = j;
printf("(%d, %d) -> %d\n", jcolor, i, j);
    }
  }

  return solver;
}


void nonlinear_solver_free(nonlinear_solver_t* solver)
{
  if (solver->initialized)
  {
    HYPRE_ParCSRBiCGSTABDestroy(solver->solver);
    HYPRE_ParCSRPilutDestroy(solver->precond);
//    HYPRE_ParCSRHybridDestroy(solver->solver);
    HYPRE_IJMatrixDestroy(solver->A);
    HYPRE_IJVectorDestroy(solver->x);
    HYPRE_IJVectorDestroy(solver->b);
    double_table_free(solver->Jij);
  }

  int num_colors = adj_graph_coloring_num_colors(solver->coloring);
  for (int color = 0; color < num_colors; ++color)
  {
    free(solver->Xps[color]);
    free(solver->Rps[color]);
    free(solver->ds[color]);
    free(solver->Jds[color]);
    free(solver->columns[color]);
  }
  free(solver->Xps);
  free(solver->Rps);
  free(solver->ds);
  free(solver->Jds);
  free(solver->columns);
  free(solver->deltaX);
  free(solver->R);
  free(solver->colors);
  adj_graph_coloring_free(solver->coloring);
  nonlinear_timestepper_free(solver->timestepper);
  nonlinear_function_free(solver->F);
  solver->index_space = NULL;
  adj_graph_free(solver->graph);
  free(solver);
}

nonlinear_function_t* nonlinear_solver_function(nonlinear_solver_t* solver)
{
  return solver->F;
}

void nonlinear_solver_compute_residual(nonlinear_solver_t* solver,
                                       double t, 
                                       double* X,
                                       double* R)
{
  ASSERT(X != NULL);
  ASSERT(R != NULL);
  nonlinear_function_eval(solver->F, t, X, R);
}

// This function estimates the product of the Jacobian (about a reference 
// state X) with the given vector y using finite differences.
static void compute_Jd_with_finite_differences(nonlinear_solver_t* solver,
                                               double t,
                                               double* X, 
                                               int color)
{
  ASSERT(X != NULL);

  int num_vertices = adj_graph_num_vertices(solver->graph);

  // Set up our work arrays for this color.
  double* R0 = solver->R;
  double* d = solver->ds[color];
  double* Jd = solver->Jds[color];
  double* Xp = solver->Xps[color];
  double* Rp = solver->Rps[color];

  // Copy the reference state to our work array.
  memcpy(Xp, X, sizeof(double) * num_vertices);

  // Perturb the the state at each of the vertices with this color.
  int pos = 0, i;
  while (adj_graph_coloring_next_vertex(solver->coloring, color, &pos, &i))
  {
    double dX = (X[i] == 0.0) ? solver->delta : solver->delta * X[i];
    Xp[i] = X[i] + dX * d[i];
  }
    
  // Compute the residual using the perturbed state.
  nonlinear_function_eval(solver->F, t, Xp, Rp);

  // Traverse the rows of the Jacobian and form the matrix project J * d.
  memset(Jd, 0, sizeof(double) * num_vertices);
  for (int i = 0; i < num_vertices; ++i)
  {
    double dX = (Xp[i] - X[i]);
    if (dX == 0.0) continue;

    // Diagonal term.
    {
      double dR = (Rp[i] - R0[i]);
      Jd[i] += dR/dX * d[i];
    }

    // Off-diagonal terms.
    int num_edges = adj_graph_num_edges(solver->graph, i);
    int *edges = adj_graph_edges(solver->graph, i);
    for (int e = 0; e < num_edges; ++e)
    {
      int j = edges[e];
      double dR = (Rp[j] - R0[j]);
printf("J(%d, %d) = %g/%g = %g\n", i, j, dR, dX, dR/dX*d[j]);
      Jd[i] += dR/dX * d[j];
    }
  }
}

void nonlinear_solver_compute_jacobian(nonlinear_solver_t* solver,
                                       double t,
                                       double* X,
                                       double_table_t* Jij)
{
  ASSERT(X != NULL);
  ASSERT(Jij != NULL);

  // Compute the residual in the reference state.
  nonlinear_function_eval(solver->F, t, X, solver->R);

  // Coleman et. al (1983) explain how a set of columns of a given color can 
  // be evaluated with a single evaluation of the residual function, so the 
  // number of evaluations of the residual function is equal to the number 
  // of colors times the number of solution components.

  // Iterate over the colors of the graph and compute the components of 
  // the Jacobian. 
  int num_colors = adj_graph_coloring_num_colors(solver->coloring);
  int num_vertices = solver->num_sites;
  for (int color = 0; color < num_colors; ++color)
  {
    // Compute the vector product J * d for this color.
    compute_Jd_with_finite_differences(solver, t, X, color);
  }

  // Fill in the entries of the Jacobian, stepping over the column
  // vectors Jd and reading off their rows Aij.
  for (int color = 0; color < num_colors; ++color)
  {
    double* d = solver->ds[color];
    double* Jd = solver->Jds[color];
    int* columns = solver->columns[color];

    // Go over the rows of this column vector and figure out the column in 
    // the Jacobian for this row.
    for (int i = 0; i < num_vertices; ++i)
    {
      int j = columns[i];
      if (j != -1)
      {
printf("column for color %d and row %d is %d\n", color, i, j);
        double Aij = Jd[i];
        printf("J(%d, %d) is %g\n", i, j, Aij);
        double_table_insert(Jij, i, j, Aij);
      }
    }
  }
#if 0
    int num_edges = adj_graph_num_edges(solver->graph, i);
    int *edges = adj_graph_edges(solver->graph, i);
    for (int e = 0; e < num_edges; ++e)
    {
      int j = edges[e];

      // What color is vertex j? 
      int color = solver->colors[j];

      // Read the (i,j)th matrix entry from the proper J*d vector.
      double Aij = solver->Jds[color][i];
printf("J(%d, %d) is %g\n", i, j, Aij);
      double_table_insert(Jij, i, j, Aij);
    }
#endif
}

static void initialize(nonlinear_solver_t* solver)
{
  if (!solver->initialized)
  {
    solver->A = HYPRE_IJMatrixNew(solver->index_space);
    solver->x = HYPRE_IJVectorNew(solver->index_space);
    solver->b = HYPRE_IJVectorNew(solver->index_space);
    // Use HYPRE's BiCGSTAB solver with the Parallel ILU (Pilut) preconditioner.
    HYPRE_ParCSRBiCGSTABCreate(adj_graph_comm(solver->graph), &solver->solver);
    HYPRE_ParCSRPilutCreate(adj_graph_comm(solver->graph), &solver->precond);
    HYPRE_ParCSRPilutSetDropTolerance(solver->precond, 1e-6);
    HYPRE_ParCSRBiCGSTABSetPrecond(solver->solver, HYPRE_ParCSRPilutSolve, 
                                   HYPRE_ParCSRPilutSetup, solver->precond);
    // Uncomment this line to see the solver's narrative.
    HYPRE_ParCSRBiCGSTABSetPrintLevel(solver->solver, 1);

#if 0
    // Use HYPRE's hybrid multigrid/Krylov solver.
    HYPRE_ParCSRHybridCreate(&solver->solver);
    // Solver types: 1 is PCG, 2 is GMRES, 3 is BiCGSTAB.
    HYPRE_ParCSRHybridSetSolverType(solver->solver, 3);
    // Set the dimension of the Krylov subspace.
    HYPRE_ParCSRHybridSetKDim(solver->solver, 10);
    // Uncomment this line to see the solver's narrative.
    HYPRE_ParCSRHybridSetPrintLevel(solver->solver, 1);
#endif

    solver->Jij = double_table_new();
    solver->initialized = true;
  }
}

static void solve_linearized_system(nonlinear_solver_t* solver, HYPRE_IJMatrix A, HYPRE_IJVector b, HYPRE_IJVector x)
{
HYPRE_IJMatrixPrint(solver->A, "A");
HYPRE_IJVectorPrint(solver->b, "b");
HYPRE_IJVectorPrint(solver->x, "x");
  HYPRE_ParCSRMatrix Aobj;
  int err = HYPRE_IJMatrixGetObject(solver->A, (void**)&Aobj);
  ASSERT(err == 0);
  ASSERT(Aobj != NULL);
  HYPRE_ParVector xobj, bobj;
  err = HYPRE_IJVectorGetObject(solver->x, (void**)&xobj);
  ASSERT(err == 0);
  ASSERT(xobj != NULL);
  err = HYPRE_IJVectorGetObject(solver->b, (void**)&bobj);
  ASSERT(err == 0);
  ASSERT(bobj != NULL);

  err = HYPRE_ParCSRBiCGSTABSetup(solver->solver, Aobj, bobj, xobj);
  ASSERT(err == 0);
  err = HYPRE_ParCSRBiCGSTABSolve(solver->solver, Aobj, bobj, xobj);
#if 0
  err = HYPRE_ParCSRHybridSetup(solver->solver, Aobj, bobj, xobj);
  ASSERT(err == 0);
  err = HYPRE_ParCSRHybridSolve(solver->solver, Aobj, bobj, xobj);
#endif
  if (err == HYPRE_ERROR_CONV)
    polymec_error("nonlinear_solver: Linear solve did not converge.");

  HYPRE_ClearAllErrors();
}

void nonlinear_solver_step(nonlinear_solver_t* solver,
                           double* t, double* x)
{
  // Make sure the solver is initialized.
  initialize(solver);

  int num_sites = (solver->index_space->high - solver->index_space->low);

  // Choose an initial timestep.
  char* explanation;
  double dt = nonlinear_timestepper_step_size(solver->timestepper, &explanation);
  log_detail("nonlinear_solver_step: Chose dt = %g (%s)", dt, explanation);

  // Compute the residual vector at the unperturbed state.
  nonlinear_solver_compute_residual(solver, *t, x, solver->R);

  bool converged = false;
  int num_iters = 0;
  double last_dt = dt;
  while (!converged && (num_iters <= solver->max_iters))
  {
    // Here we go.
    ++num_iters;

    // Choose a new timestep if needed.
    dt = nonlinear_timestepper_step_size(solver->timestepper, &explanation);
    if (last_dt != dt)
    {
      log_detail("nonlinear_solver_step: Chose dt = %g (%s)", dt, explanation);
      last_dt = dt;
    }

    // Compute the components of the Jacobian as needed.
    if (nonlinear_timestepper_recompute_jacobian(solver->timestepper))
    {
      log_detail("nonlinear_solver_step: Computing Jacobian at t = %g.", *t);
      nonlinear_solver_compute_jacobian(solver, *t, x, solver->Jij);
    }

    // Compute J = I/dt - dR/dX to implement the backward Euler method.
    double_table_t* backward_euler_J = double_table_new();
    double_table_cell_pos_t pos = double_table_start(solver->Jij);
    int i, j;
    double Jij;
    while (double_table_next_cell(solver->Jij, &pos, &i, &j, &Jij))
    {
      if (i == j)
        double_table_insert(backward_euler_J, i, j, 1.0/dt - Jij);
      else
        double_table_insert(backward_euler_J, i, j, -Jij);
    }

    // Set up the linear system.
    HYPRE_IJMatrixSetValuesFromTable(solver->A, solver->index_space, backward_euler_J);
    HYPRE_IJVectorSetValuesFromArray(solver->b, solver->index_space, solver->R);
    double_table_free(backward_euler_J);

    // Now solve the linear system for the increment.
    solve_linearized_system(solver, solver->A, solver->b, solver->x);

    // Copy the solution to dX.
    HYPRE_IJVectorGetValuesToArray(solver->x, solver->index_space, solver->deltaX);

    // Add the increment to x.
    for (int i = 0; i < num_sites; ++i)
      x[i] += solver->deltaX[i];

    // Recompute the residual and measure its norm.
    nonlinear_solver_compute_residual(solver, *t, x, solver->R);
    double res_norm = 0.0;
    for (int i = 0; i < num_sites; ++i)
    {
      double Ri = solver->R[i];
      res_norm += Ri*Ri;
    }
    res_norm = sqrt(res_norm);

    // Is the norm small enough?
    if (res_norm < solver->tolerance)
    {
      // Success!
      nonlinear_timestepper_converged(solver->timestepper, dt, 
                                      res_norm, num_iters);
      converged = true;
    }
    else
    {
      // Awwwww...
      nonlinear_timestepper_failed(solver->timestepper, dt);
    }
  }

  // If we exceeded the maximum number of iterations, we're toast.
  if (num_iters > solver->max_iters)
  {
    polymec_error("nonlinear_solver_step: Maximum number of iterations (%d) exceeded", 
                  solver->max_iters);
  }

  // Update the time.
  *t += dt;
}

void nonlinear_solver_solve(nonlinear_solver_t* solver, double t, double* x)
{
  ASSERT(solver->timestepper == NULL);

  // Make sure the solver is initialized.
  initialize(solver);

  int num_sites = (solver->index_space->high - solver->index_space->low);

  // Compute the residual vector at the unperturbed state.
  nonlinear_solver_compute_residual(solver, t, x, solver->R);

  bool converged = false;
  int num_iters = 0;
  while (!converged && (num_iters <= solver->max_iters))
  {
    // Here we go.
    ++num_iters;

    // Compute the components of the Jacobian.
    nonlinear_solver_compute_jacobian(solver, t, x, solver->Jij);

    // Set up the linear system.
    HYPRE_IJMatrixSetValuesFromTable(solver->A, solver->index_space, solver->Jij);
    HYPRE_IJVectorSetValuesFromArray(solver->b, solver->index_space, solver->R);

    // Now solve the linear system for the increment.
    solve_linearized_system(solver, solver->A, solver->b, solver->x);

    // Copy the solution to deltaX.
    HYPRE_IJVectorGetValuesToArray(solver->x, solver->index_space, solver->deltaX);

    // Add the increment to x.
    for (int i = 0; i < num_sites; ++i)
      x[i] += solver->deltaX[i];

    // Recompute the residual and measure its norm.
    nonlinear_solver_compute_residual(solver, t, x, solver->R);
    double res_norm = 0.0;
    for (int i = 0; i < num_sites; ++i)
    {
      double Ri = solver->R[i];
      res_norm += Ri*Ri;
    }
    res_norm = sqrt(res_norm);

    // Is the norm small enough?
    if (res_norm < solver->tolerance)
      converged = true;
  }

  // If we exceeded the maximum number of iterations, we're toast.
  if (num_iters > solver->max_iters)
  {
    polymec_error("nonlinear_solver_solve: Maximum number of iterations (%d) exceeded", 
                  solver->max_iters);
  }
}

void nonlinear_solver_integrate(nonlinear_solver_t* solver,
                                double t1, double* x1,
                                double t2, double* x2)
{
  ASSERT(solver->timestepper != NULL);
  ASSERT(x1 != NULL);
  ASSERT(x2 != NULL);

  double t = t1;
  int N = solver->index_space->high - solver->index_space->low;
  memcpy(x2, x1, sizeof(double) * N);
  while (t < t2)
    nonlinear_solver_step(solver, &t, x2);
}

double nonlinear_solver_tolerance(nonlinear_solver_t* solver)
{
  return solver->tolerance;
}
                                  
void nonlinear_solver_set_tolerance(nonlinear_solver_t* solver, double tolerance)
{
  ASSERT(tolerance > 0.0);
  solver->tolerance = tolerance;
}

int nonlinear_solver_max_iterations(nonlinear_solver_t* solver)
{
  return solver->max_iters;
}
                                  
void nonlinear_solver_set_max_iterations(nonlinear_solver_t* solver, int max_iters)
{
  ASSERT(max_iters > 1);
  solver->max_iters = max_iters;
}
                                  
