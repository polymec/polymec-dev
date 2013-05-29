#include "integrators/nonlinear_solver.h"

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
  double* XX; // Work array
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
  solver->num_sites = graph->num_vertices;
  int N = nonlinear_function_num_comps(F);
  solver->XX = malloc(sizeof(double) * solver->num_sites * N);
  return solver;
}


void nonlinear_solver_free(nonlinear_solver_t* solver)
{
  nonlinear_function_free(solver->F);
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
  // FIXME
}

