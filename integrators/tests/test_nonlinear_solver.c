#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmockery.h"
#include "core/st_func.h"
#include "integrators/nonlinear_solver.h"
#include "integrators/simple_nonlinear_timestepper.h"

// This function constructs a block adjacency graph for a set of cells in 
// 1D.
static adj_graph_t* block_graph_1d(int num_cells, int block_size)
{
  // The adjacency graph for a 1D grid is fairly simple.
  adj_graph_t* graph = adj_graph_new(MPI_COMM_WORLD, num_cells);
  for (int v = 0; v < adj_graph_num_vertices(graph); ++v)
  {
    if (v == 0)
    {
      adj_graph_set_num_edges(graph, v, 1);
      int* e = adj_graph_edges(graph, v);
      e[0] = 1;
    }
    else if (v == (num_cells-1))
    {
      adj_graph_set_num_edges(graph, v, 1);
      int* e = adj_graph_edges(graph, v);
      e[0] = num_cells-2;
    }
    else
    {
      adj_graph_set_num_edges(graph, v, 2);
      int* e = adj_graph_edges(graph, v);
      e[0] = v - 1;
      e[1] = v + 1;
    }
  }

  if (block_size > 1)
  {
    // Create an adjacency graph that reflects the connectivity of 
    // components of the solution.
    adj_graph_t* block_graph = adj_graph_new_with_block_size(block_size, graph);
    adj_graph_free(graph);
    return block_graph;
  }
  else 
    return graph;
}

// This function populates the given vector X with the values of the 
// function at the time t, evaluated at N equally-spaced cells spanning 
// the interval [x1, x2].
static void initialize_solution_1d(st_func_t* F, 
                                   int N,
                                   double x1, 
                                   double x2, 
                                   double t, 
                                   double* X)
{
  ASSERT(N > 0);
  ASSERT(x1 < x2);

  double dx = (x2 - x1) / N;
  point_t x = {.x = 0.0, .y = 0.0, .z = 0.0};
  int num_comps = st_func_num_comp(F);
  for (int i = 0; i < N; ++i)
  {
    x.x = 0.5 + i * dx;
    st_func_eval(F, &x, t, &X[num_comps*i]);
  }
}

// This function maps a solution vector X to itself.
static void identity_transformation(void* context,
                                    double t, 
                                    double* x,
                                    double* R)
{
  // Get the dimension.
  int N = *((int*)context); 
  
  // Do the mapping.
  memcpy(R, x, sizeof(double) * N);
}

static void test_identity_jacobian(void** state)
{
  // Set up our identity mapping.
  int N = 32;
  nonlinear_function_vtable vtable = {.eval_residual = identity_transformation};
  nonlinear_function_t* F = nonlinear_function_new("Identity mapping", 
                                                   &N, vtable); 

  // Set up a nonlinear solver.
  adj_graph_t* graph = block_graph_1d(N, 1);
  nonlinear_solver_t* solver = nonlinear_solver_new(F, NULL, graph);

  // Initialize an X vector.
  double* X = malloc(sizeof(double) * N);
  for (int i = 0; i < N; ++i)
    X[i] = 1.0 * i;

  // Compute the components of the Jacobian.
  double_table_t* Jij = double_table_new();
  nonlinear_solver_compute_jacobian(solver, 0.0, X, Jij);

  // Print the elements of the Jacobian.
  {
    double_table_cell_pos_t pos = double_table_start(Jij);
    int i, j;
    double ij;
    while (double_table_next_cell(Jij, &pos, &i, &j, &ij))
      printf("%d, %d: %g\n", i, j, ij);
  }

  // Clean up.
  double_table_free(Jij);
  nonlinear_function_free(F);
}

//------------------------------------------------------------------------
//                          1D linear advection
//------------------------------------------------------------------------
// We test our nonlinear solver on the 1D linear advection problem. 
typedef struct 
{
  int N;          // Number of grid cells.
  int num_comps;  // Number of solution components.
  double c;       // Advection velocity.
} advection_t;

// Here's the residual function for the Euler equations.
static void advection_residual(void* context,
                               double t,
                               double* x,
                               double* R)
{
  advection_t* adv = context;
  int N = adv->N;
  double dx = 1.0 / N;
  int num_comps = adv->num_comps;
  double c = adv->c;

  // Compute the inter-cell fluxes using simple upwinding.
  double fluxes[num_comps*(N+1)];
  for (int i = 1; i < N; ++i)
  {
    for (int comp = 0; comp < num_comps; ++comp)
    {
      double xL = x[num_comps*(i-1) + comp];
      double xR = x[num_comps*i + comp];

      // Upwinded solution.
      double xM;
      if (c >= 0.0) // right-moving wave
        xM = xL;
      else // left-moving wave
        xM = xR;

      // Now compute the inter-face flux.
      fluxes[num_comps*i+comp] = c * xM;
    }
  }

  // Boundary fluxes are periodic.
  for (int comp = 0; comp < num_comps; ++comp)
  {
    // Left boundary gets rightmost flux.
    fluxes[comp] = fluxes[(N-1) * num_comps + comp];
    // Right boundary gets leftmost flux.
    fluxes[N*num_comps + comp] = fluxes[num_comps + comp];
  }

  // Compute the time derivatives of the conserved quantities using 
  // the discrete divergence theorem.
  for (int i = 0; i < N; ++i)
  {
    for (int comp = 0; comp < num_comps; ++comp)
      R[num_comps*i + comp] = (fluxes[num_comps*i + comp] - fluxes[num_comps*(i+1)+comp]) / dx;
  }
}

static void advection_dtor(void* context)
{
  advection_t* adv = context;
  free(adv);
}

// Here's the method for constructing the Sod simulation.
static nonlinear_function_t* advection_new(double c, 
                                           int num_comps,
                                           int num_cells)
{
  ASSERT(num_cells > 0);

  advection_t* adv = malloc(sizeof(advection_t));
  adv->c = c;
  adv->num_comps = num_comps;
  adv->N = num_cells;

  nonlinear_function_vtable vtable = {.eval_residual = advection_residual, 
                                      .dtor = advection_dtor};
  return nonlinear_function_new("Linear advection", adv, vtable); 
}

// Here's a 1D pulse function on a periodic domain [0, 1].
static void pulse_1d(void* context,
                     point_t* x, 
                     double t,
                     double* u)
{
  advection_t* adv = context;
  double c = adv->c;
  int num_comps = adv->num_comps;

  // Find the time-shifted coordinate along the domain.
  double z = x->x - c * t;
  while (z > 1.0)
    z -= 1.0;
  while (z < 0.0)
    z += 1.0;

  // Compute the pulse height.
  for (int i = 0; i < num_comps; ++i)
    u[i] = ((z >= 0.25) && (z <= 0.75)) ? 1.0 : 0.0;
}

static void test_linear_advection(double c, 
                                  int num_comps, 
                                  int num_cells,
                                  const char* problem_name,
                                  void (*initial_cond)(void*, point_t*, double, double*),
                                  double t1,
                                  double t2)
{
  // Set up the residual function.
  nonlinear_function_t* advection_F = advection_new(c, num_comps, num_cells);

  // Set up a simple time stepper.
  double dt0 = 1.0 / num_cells;
  nonlinear_timestepper_t* stepper = simple_nonlinear_timestepper_new(dt0, 10.0*dt0, 0.5, 8, 1.2);

  // Set up the 1D adjacency graph.
  adj_graph_t* graph = block_graph_1d(num_cells, 1);

  // Create our nonlinear solver.
  nonlinear_solver_t* solver = nonlinear_solver_new(advection_F, stepper, graph);

  // Set up initial conditions.
  st_vtable X0_vtable = {.eval = initial_cond};
  st_func_t* X0 = st_func_new(problem_name, 
                              nonlinear_function_context(advection_F), 
                              X0_vtable, 
                              ST_INHOMOGENEOUS, 
                              ST_NONCONSTANT,
                              num_comps);
  double* X1 = malloc(sizeof(double) * 3 * num_cells);
  initialize_solution_1d(X0, num_cells, 0.0, 1.0, t1, X1);

  // Integrate!
  double* X2 = malloc(sizeof(double) * 3 * num_cells);
  nonlinear_solver_integrate(solver, t1, X1, t2, X2);

  // Write out the components of the solution vector X2.
  for (int i = 0; i < num_cells; ++i)
    printf("%d %g %g %g\n", i, X2[3*i], X2[3*i+1], X2[3*i+2]);

  // Clean up.
  free(X2);
  free(X1);
  nonlinear_solver_free(solver);
}

static void test_single_comp_linear_advection(void** state)
{
  test_linear_advection(1.0, 1, 32, "1D pulse", pulse_1d, 0.0, 1.0);
}

//------------------------------------------------------------------------
//                          1D Sod shock tube
//------------------------------------------------------------------------
// We test our nonlinear solver on the 1D Sod shock tube problem. The 
// following are bits of a simple 1D Eulerian hydrodynamics scheme on a 
// regular grid.
typedef struct 
{
  int N; // Number of grid cells.
  double gamma; // Polytropic index.
} sod_t;

// Here's the solution to the Riemann problem for a polytropic gas.
static void solve_riemann(double gamma,
                          double rhoL, double uL, double pL,
                          double rhoR, double uR, double pR,
                          double* rho, double* u, double* p)
{
  // Some thresholds.
  static const double small_rho = 1e-7;
  static const double small_p = 1e-7;
  static const double epsilon = 1e-8;

  // Sound speeds.
  double cL = sqrt(gamma*pL/rhoL);
  double cR = sqrt(gamma*pR/rhoR);

  // Find the pressure and velocity in the "star" region.
  double wL = rhoL * cL;
  double wR = rhoR * cR;
  double p_star = (wR*pL + wL*pR + wL*wR*(uL - uR)) / (wL + wR);
  double u_star = (wL*uL + wR*uR + pL - pR) / (wL + wR);

  // Determine which way the wave is moving and use it to reconstruct 
  // a linearized signal.
  double rho0, u0, p0, c0, sign;
  if (u_star > 0.0) // wave moving to the right
  {
    rho0 = rhoL;
    p0 = pL;
    u0 = uL;
    c0 = cL;
    sign = 1.0;
  }
  else // wave moving to the left
  {
    rho0 = rhoR;
    p0 = pR;
    u0 = uR;
    c0 = cR;
    sign = -1.0;
  }

  // Determine more quantities in the star region.
  double rho_star = MAX(rho0 + (p_star - p0) / (c0*c0), small_rho);
  double c_star = sqrt(fabs(gamma * p_star/rho_star));
  double w_star = 0.5 * (c_star * rho_star + c0 * rho0);

  // Construct the linearized signal.
  double sp_out = c0 - sign * u0;
  double sp_in = c_star - sign * u_star;
  double u_shock = w_star/rho_star - sign * u_star;
  if (p_star > p0)
  {
    sp_out = u_shock;
    sp_in = u_shock;
  }

  double frac = ((1.0 + (sp_out + sp_in)/MAX(sp_out - sp_in, epsilon))/2.0);
  frac = MAX(0.0, MIN(1.0, frac));
  *rho = rho0 + frac * (rho_star - rho0);
  *u = u0 + frac * (u_star - u0);
  *p = p0 + frac * (p_star - p0);

  // Apply limits.
  if (sp_out <= 0.0)
  {
    *rho = rho0;
    *u = u0;
    *p = p0;
  }
  if (sp_in >= 0.0)
  {
    *rho = rho_star;
    *u = u_star;
    *p = p_star;
  }
  *rho = MAX(*rho, small_rho);
  *p = MAX(*p, small_p);

}

// Here's the residual function for the Euler equations.
static void sod_residual(void* context,
                         double t,
                         double* x,
                         double* R)
{
  sod_t* sod = context;
  int N = sod->N;
  double gamma = sod->gamma;
  double dx = 1.0 / N;

  // Compute the inter-cell fluxes by solving Riemann problems.
  double fluxes[3*(N+1)];
  for (int i = 1; i < N; ++i)
  {
    // Left-side variables.
    double rhoL = x[3*(i-1)];
    double uL   = x[3*(i-1)+1] / rhoL;
    double EL   = x[3*(i-1)+2] / rhoL;
    double pL = (gamma - 1.0) * rhoL * (EL - 0.5 * uL * uL);

    double rhoR = x[3*i];
    double uR   = x[3*i+1] / rhoR;
    double ER   = x[3*i+2] / rhoR;
    double pR = (gamma - 1.0) * rhoR * (ER - 0.5 * uR * uR);

    // Solve the Riemann problem at the left and right faces to obtain
    // values for the fluxes.
    double rho, u, p;
    solve_riemann(gamma, rhoL, uL, pL, rhoR, uR, pR, &rho, &u, &p);

    // Now compute the inter-face flux.
    fluxes[3*i]   = rho*u;
    fluxes[3*i+1] = rho*u*u + p;
    fluxes[3*i+2] = rho*u*(0.5*u*u + p);
  }

  // Boundary fluxes are zero, except pressures, which are kept constant.
  fluxes[0] = fluxes[2] = 0.0;
  double rhoE0 = x[2];
  double p0 = (gamma - 1.0) * rhoE0;
  fluxes[1] = p0;

  fluxes[3*N] = fluxes[3*N+2] = 0.0;
  double rhoEN = x[3*(N-1)+2];
  double pN = (gamma - 1.0) * rhoEN;
  fluxes[3*N+1] = pN;

  // Compute the time derivatives of the conserved quantities using 
  // the discrete divergence theorem.
  for (int i = 0; i < N; ++i)
  {
    R[3*i]   = (fluxes[3*i] - fluxes[3*(i+1)])/dx;
    R[3*i+1] = (fluxes[3*i+1] - fluxes[3*(i+1)+1])/dx;
    R[3*i+2] = (fluxes[3*i+2] - fluxes[3*(i+1)+2])/dx;
  }
}

static void sod_dtor(void* context)
{
  sod_t* sod = context;
  free(sod);
}

// Here's the method for constructing the Sod simulation.
static nonlinear_function_t* sod_new(double gamma, 
                                     int num_cells)
{
  ASSERT(gamma > 1.0);
  ASSERT(num_cells > 0);

  sod_t* sod = malloc(sizeof(sod_t));
  sod->gamma = gamma;
  sod->N = num_cells;

  nonlinear_function_vtable vtable = {.eval_residual = sod_residual, 
                                      .dtor = sod_dtor};
  return nonlinear_function_new("Sod shock tube", sod, vtable); 
}

static void test_sod_shock_tube(void** state)
{
  double gamma = 5.0/3.0;
  double rho1 = 1.0, p1 = 1.0;
  double rho2 = 0.125, p2 = 0.1;
  double t_final = 0.2;
  int N = 32;
  nonlinear_function_t* sod_F = sod_new(gamma, N);

  double dt0 = 1.0/N;
  nonlinear_timestepper_t* stepper = simple_nonlinear_timestepper_new(dt0, 10.0*dt0, 0.5, 8, 1.2);

  // Set up the 1D adjacency graph.
  adj_graph_t* graph = block_graph_1d(N, 3);

  // Create a solution array (rho, rho*u, rho*E) and initialize it.
  double dx = 1.0 / N;
  double* X1 = malloc(sizeof(double) * 3 * N);
  double rhoE1 = p1 / (gamma - 1.0);
  double rhoE2 = p2 / (gamma - 1.0);
  for (int i = 0; i < N; ++i)
  {
    double x = (0.5 + i) * dx;
    if (x < 0.5)
    {
      X1[3*i]   = rho1;
      X1[3*i+2] = rhoE1;
    }
    else
    {
      X1[3*i]   = rho2;
      X1[3*i+2] = rhoE2;
    }
    X1[3*i+1] = 0.0;
  }

  // Create our nonlinear solver.
  nonlinear_solver_t* solver = nonlinear_solver_new(sod_F, stepper, graph);

  // Integrate!
  double* X2 = malloc(sizeof(double) * 3 * N);
  nonlinear_solver_integrate(solver, 0.0, X1, t_final, X2);

  // Write out the components of the solution vector X2.
  for (int i = 0; i < N; ++i)
    printf("%d %g %g %g\n", i, X2[3*i], X2[3*i+1], X2[3*i+2]);

  // Clean up.
  free(X2);
  free(X1);
  nonlinear_solver_free(solver);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  set_log_level(LOG_DEBUG);
  const UnitTest tests[] = 
  {
    unit_test(test_identity_jacobian)
//    unit_test(test_single_comp_linear_advection),
//    unit_test(test_sod_shock_tube)
  };
  return run_tests(tests);
}
