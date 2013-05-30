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

struct nonlinear_solver_t 
{
  nonlinear_function_t* F;
  adj_graph_t* graph;
  int num_sites;
  double* XX; // Work array.

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

// Returns the context object for this nonlinear function.
void* nonlinear_function_context(nonlinear_function_t* F)
{
  return F->context;
}

// Returns the number of components per site in this nonlinear function.
int nonlinear_function_num_comps(nonlinear_function_t* F)
{
  return F->num_comps;
}

void nonlinear_function_eval(nonlinear_function_t* F,
                             int site,
                             double* X,
                             double* R)
{
  ASSERT(site >= 0);
  ASSERT(X != NULL);
  ASSERT(R != NULL);
  F->vtable.eval_residual(F->context, F->num_comps, site, X, R);
}

nonlinear_solver_t* nonlinear_solver_new(nonlinear_function_t* F,
                                         adj_graph_t* graph)
{
  ASSERT(F != NULL);
  ASSERT(num_sites > 0);

  nonlinear_solver_t* solver = malloc(sizeof(nonlinear_solver_t));
  solver->F = F;
  solver->graph = graph;
  solver->num_sites = adj_graph_num_vertices(graph);
  int N = nonlinear_function_num_comps(F);
  solver->XX = malloc(sizeof(double) * solver->num_sites * N);
  solver->index_space = index_space_new_from_low_and_high(adj_graph_comm(graph),
                                                          adj_graph_first_vertex(graph),
                                                          adj_graph_last_vertex(graph));
  solver->initialized = false;
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
  }

  nonlinear_function_free(solver->F);
  solver->index_space = NULL;
  adj_graph_free(solver->graph);
  free(solver->XX);
  free(solver);
}

nonlinear_function_t* nonlinear_solver_function(nonlinear_solver_t* solver)
{
  return solver->F;
}

void nonlinear_solver_compute_residual(nonlinear_solver_t* solver,
                                       int site,
                                       double* X,
                                       double* R)
{
  ASSERT(site >= 0);
  ASSERT(site < solver->num_sites);
  ASSERT(X != NULL);
  ASSERT(R != NULL);
  nonlinear_function_eval(solver->F, site, X, R);
}

void nonlinear_solver_compute_jacobian(nonlinear_solver_t* solver,
                                       int site,
                                       double* X,
                                       double delta,
                                       double_table_t* Jij)
{
  ASSERT(site >= 0);
  ASSERT(site < solver->num_sites);
  ASSERT(X != NULL);
  ASSERT(delta > 0.0);
  ASSERT(Jij != NULL);

  // Compute the residual in the reference state.
  int N = nonlinear_function_num_comps(solver->F);
  double R0[N];
  nonlinear_function_eval(solver->F, site, X, R0);

  // Increment the various solution coefficients to compute the components 
  // of the Jacobian. We only need to increment those coefficients that 
  // are connected to this site in the adjacency graph.
  for (int i = 0; i < N; ++i) 
  {
    // Copy the reference state to our work array.
    memcpy(solver->XX, X, sizeof(double) * N * solver->num_sites);

    // Find the other sites that affect this one, and build a list of all 
    // sites to be incremented.
    int num_adj_sites = adj_graph_num_edges(solver->graph, site);
    int inc_sites[1+num_adj_sites], k = 1;
    inc_sites[0] = site;
    int num_edges = adj_graph_num_edges(solver->graph, site);
    int* edges = adj_graph_edges(solver->graph, site);
    for (int j = 0; j < num_edges; ++j)
      inc_sites[k] = edges[j];

    // Increment the ith component of each of these sites within the work array.
    double dXs[1 + num_adj_sites];
    for (int j = 0; j < 1 + num_adj_sites; ++j)
    {
      int k = inc_sites[j];
      dXs[j] = (X[k] == 0.0) ? delta : delta * X[k];
      solver->XX[k] = X[k] + dXs[j];
    }

    // Compute the residual for the incremented solution vector.
    double R[N];
    nonlinear_function_eval(solver->F, site, solver->XX, R);

    // Stash the resulting partial derivatives in Jij.
    for (int j = 0; j < N; ++j)
    {
      for (int k = 0; k < 1 + num_adj_sites; ++k)
      {
        int l = inc_sites[k];
        //       dR    (R - R0)
        // Jij = --j = -------- , 
        //       dX       dX
        //         i
        //
        // and remember that we are computing the components of the 
        // entire Jacobian, so the "i" up here is the index of the 
        // ith component at the incremented site, N*l+i.
        // The "j" is the jth component of the given site in the residual, 
        // or N*site + j.
        double_table_insert(Jij, N*l+i, N*site+j, (R[j] - R0[j])/dXs[k]);
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
    solver->initialized = true;
  }
}

static inline void solve(nonlinear_solver_t* solver, HYPRE_IJMatrix A, HYPRE_IJVector b, HYPRE_IJVector x)
{
//HYPRE_IJMatrixPrint(solver->A, "A");
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
    polymec_error("nonlinear_solver: Solve did not converge.");

  HYPRE_ClearAllErrors();
}

void nonlinear_solver_step(nonlinear_solver_t* solver,
                           double* t, double* x)
{
  // Make sure the solver is initialized.
  initialize(solver);

  static double epsilon = FLT_EPSILON * FLT_MIN; // Smallest possible float.

  // Compute the components of the Jacobian and the residual at each site.
  int N = solver->index_space->high - solver->index_space->low;
  int num_comps = nonlinear_function_num_comps(solver->F);
  for (int i = 0; i < N; ++i)
  {
    nonlinear_solver_compute_jacobian(solver, i, x, epsilon, solver->Jij);
    nonlinear_solver_compute_residual(solver, i, x, solver->R);
  }

  // Set up the linear system.
  HYPRE_IJMatrixSetValuesFromTable(solver->A, solver->index_space, solver->Jij);
//HYPRE_IJMatrixPrint(solver->A, "L");
  HYPRE_IJVectorSetValuesFromArray(solver->b, solver->index_space, solver->R);

  // Now solve the linear system for the increment.
  solve(solver, solver->A, solver->b, solver->x);

  // Copy the solution to XX.
  HYPRE_IJVectorGetValuesToArray(solver->x, solver->index_space, solver->XX);

  // Add the increment to x.
  for (int i = 0; i < N*num_comps; ++i)
    x[i] += solver->XX[i];
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

