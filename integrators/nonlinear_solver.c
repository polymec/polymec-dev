#include <float.h>
#include "integrators/nonlinear_solver.h"
#include "core/hypre_helpers.h"

struct nonlinear_function_t 
{
  void* context;
  char* name;
  int num_comps;
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
  adj_graph_t* graph;
  adj_graph_coloring_t* coloring;
  int num_sites;
  double *XX, *RR; // Work arrays.

  // Linear system and solver.
  HYPRE_IJMatrix A;
  HYPRE_IJVector x, b;
  HYPRE_Solver solver;

  // Flag is set to true if the linear system above is initialized.
  bool initialized;

  // Index space.
  index_space_t* index_space;

  // Jacobian and residual coefficients.
  double_table_t* Jij;
  double* R;

  // Residual norm tolerance.
  double tolerance;

  // Maximum number of iterations before abject failure.
  int max_iters;
};

nonlinear_function_t* nonlinear_function_new(const char* name, 
                                             void* context,
                                             int num_comps,
                                             nonlinear_function_vtable vtable)
{
  ASSERT(num_comps >= 1);
  ASSERT(vtable.eval_residual != NULL);

  nonlinear_function_t* F = malloc(sizeof(nonlinear_function_t));
  F->name = strdup(name);
  F->context = context;
  F->num_comps = num_comps;
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

int nonlinear_function_num_comps(nonlinear_function_t* F)
{
  return F->num_comps;
}

void nonlinear_function_eval(nonlinear_function_t* F,
                             double t,
                             double* X,
                             double* R)
{
  ASSERT(X != NULL);
  ASSERT(R != NULL);
  F->vtable.eval_residual(F->context, F->num_comps, t, X, R);
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
  ASSERT(timestepper != NULL);
  ASSERT(graph != NULL);

  nonlinear_solver_t* solver = malloc(sizeof(nonlinear_solver_t));
  solver->F = F;
  solver->timestepper = timestepper;
  solver->graph = graph;
  solver->coloring = NULL;
  solver->num_sites = adj_graph_num_vertices(graph);
  int num_comps = nonlinear_function_num_comps(F);
  solver->XX = malloc(sizeof(double) * solver->num_sites * num_comps);
  solver->RR = malloc(sizeof(double) * solver->num_sites * num_comps);
  solver->index_space = index_space_from_low_and_high(adj_graph_comm(graph),
                                                      adj_graph_first_vertex(graph),
                                                      adj_graph_last_vertex(graph)+1);
  solver->initialized = false;
  solver->tolerance = 1e-8;
  solver->max_iters = 100;
  return solver;
}


void nonlinear_solver_free(nonlinear_solver_t* solver)
{
  if (solver->initialized)
  {
    HYPRE_ParCSRHybridDestroy(solver->solver);
    HYPRE_IJMatrixDestroy(solver->A);
    HYPRE_IJVectorDestroy(solver->x);
    HYPRE_IJVectorDestroy(solver->b);
    double_table_free(solver->Jij);
    free(solver->R);
    adj_graph_coloring_free(solver->coloring);
  }

  nonlinear_timestepper_free(solver->timestepper);
  nonlinear_function_free(solver->F);
  solver->index_space = NULL;
  adj_graph_free(solver->graph);
  free(solver->RR);
  free(solver->XX);
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

void nonlinear_solver_compute_jacobian(nonlinear_solver_t* solver,
                                       double t,
                                       double* X,
                                       double delta,
                                       double_table_t* Jij)
{
  ASSERT(X != NULL);
  ASSERT(delta > 0.0);
  ASSERT(Jij != NULL);

  // Compute the residual in the reference state.
  int num_comps = nonlinear_function_num_comps(solver->F);
  nonlinear_function_eval(solver->F, t, X, solver->R);

  // Iterate over the colors of the graph and compute the components of 
  // the Jacobian.
  int num_colors = adj_graph_coloring_num_colors(solver->coloring);
  int num_sites = solver->num_sites;
  for (int color = 0; color < num_colors; ++color)
  {
    // Increment the various solution coefficients to compute the components 
    // of the Jacobian. 
    for (int inc_comp = 0; inc_comp < num_comps; ++inc_comp) 
    {
      // Copy the reference state to our work array.
      memcpy(solver->XX, X, sizeof(double) * num_comps * solver->num_sites);

      // Increment this component at each of the sites.
      for (int site = 0; site < num_sites; ++site)
      {
        int k = num_comps*site + inc_comp;
        double dX = (X[k] == 0.0) ? delta : delta * X[k];
        solver->XX[k] = X[k] + dX;
      }

      // Compute the residual for the incremented solution vector.
      nonlinear_function_eval(solver->F, t, solver->XX, solver->RR);

      // Stash the resulting partial derivatives for each of the vertices 
      // belonging to this color.
      int pos = 0, vertex;
printf("color: %d\n", color);
      while (adj_graph_coloring_next_vertex(solver->coloring, color, &pos, &vertex))
      {
printf("  vertex: %d\n", vertex);
        for (int res_comp = 0; res_comp < num_comps; ++res_comp)
        {
          // Column of Jacobian = component of residual vector.
          int col = num_comps * vertex + res_comp;
          for (int site = 0; site < num_sites; ++site)
          {
            // Row of Jacobian = component of solution vector.
            int row = num_comps * site + inc_comp;

            //       dR    (R - R0)
            // Jij = --j = -------- , 
            //       dX       dX
            //         i
            double dX = solver->XX[row] - X[row];
            double_table_insert(Jij, row, col, (solver->RR[col] - solver->R[col])/dX);
          }
        }
      }
    }
  }
}

static void initialize(nonlinear_solver_t* solver)
{
  if (!solver->initialized)
  {
    solver->A = HYPRE_IJMatrixNew(solver->index_space);
    solver->x = HYPRE_IJVectorNew(solver->index_space);
    solver->b = HYPRE_IJVectorNew(solver->index_space);
    HYPRE_ParCSRHybridCreate(&solver->solver);
    HYPRE_ParCSRHybridSetSolverType(solver->solver, 2);
    HYPRE_ParCSRHybridSetKDim(solver->solver, 5);
    solver->Jij = double_table_new();
    int N = nonlinear_function_num_comps(solver->F);
    solver->R = malloc(sizeof(double) * N * adj_graph_num_vertices(solver->graph));
    solver->coloring = adj_graph_coloring_new(solver->graph, 
                                              SMALLEST_LAST); 
    log_detail("nonlinear_solver: graph coloring produced %d colors.", 
               adj_graph_coloring_num_colors(solver->coloring));
    solver->initialized = true;
  }
}

static inline void solve_linearized_system(nonlinear_solver_t* solver, HYPRE_IJMatrix A, HYPRE_IJVector b, HYPRE_IJVector x)
{
HYPRE_IJMatrixPrint(solver->A, "A");
//HYPRE_IJVectorPrint(solver->b, "b");
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

  err = HYPRE_ParCSRHybridSetup(solver->solver, Aobj, bobj, xobj);
  ASSERT(err == 0);
  err = HYPRE_ParCSRHybridSolve(solver->solver, Aobj, bobj, xobj);
  if (err == HYPRE_ERROR_CONV)
    polymec_error("nonlinear_solver: Linear solve did not converge.");

  HYPRE_ClearAllErrors();
}

void nonlinear_solver_step(nonlinear_solver_t* solver,
                           double* t, double* x)
{
  // Make sure the solver is initialized.
  initialize(solver);

  static double epsilon = FLT_EPSILON; //  * FLT_MIN; // Smallest possible float.

  int num_sites = solver->index_space->high - solver->index_space->low;
  int num_comps = nonlinear_function_num_comps(solver->F);

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
      nonlinear_solver_compute_jacobian(solver, *t, x, epsilon, solver->Jij);
    }

    // Go over the values in the table and compute 1.0/dt - Jij to 
    // implement the backward Euler method.
    double_table_t* be_J = double_table_new();
    double_table_cell_pos_t pos = double_table_start(solver->Jij);
    int i, j;
    double Jij;
    while (double_table_next_cell(solver->Jij, &pos, &i, &j, &Jij))
      double_table_insert(be_J, i, j, 1.0/dt - Jij);

    // Set up the linear system.
    HYPRE_IJMatrixSetValuesFromTable(solver->A, solver->index_space, be_J);
    //HYPRE_IJMatrixPrint(solver->A, "L");
    HYPRE_IJVectorSetValuesFromArray(solver->b, solver->index_space, solver->R);
    double_table_free(be_J);

    // Now solve the linear system for the increment.
    solve_linearized_system(solver, solver->A, solver->b, solver->x);

    // Copy the solution to XX.
    HYPRE_IJVectorGetValuesToArray(solver->x, solver->index_space, solver->XX);

    // Add the increment to x.
    for (int i = 0; i < num_sites*num_comps; ++i)
      x[i] += solver->XX[i];

    // Recompute the residual and measure its norm.
    nonlinear_solver_compute_residual(solver, *t, x, solver->R);
    double res_norm = 0.0;
    for (int i = 0; i < num_comps * num_sites; ++i)
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

void nonlinear_solver_integrate(nonlinear_solver_t* solver,
                                double t1, double* x1,
                                double t2, double* x2)
{
  double t = t1;
  int N = solver->index_space->high - solver->index_space->low;
  int num_comps = nonlinear_function_num_comps(solver->F);
  memcpy(x2, x1, sizeof(double) * N * num_comps);
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
                                  
