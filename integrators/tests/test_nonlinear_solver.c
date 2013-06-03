#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmockery.h"
#include "integrators/nonlinear_solver.h"
#include "integrators/simple_nonlinear_timestepper.h"

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
                         int num_comps,
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
  for (int i = 1; i < N-1; ++i)
  {
    // Left-side variables.
    double rhoL = x[num_comps*i];
    double uL   = x[num_comps*i+1] / rhoL;
    double EL   = x[num_comps*i+2] / rhoL;
    double pL = (gamma - 1.0) * rhoL * (EL - 0.5 * rhoL * uL * uL);

    double rhoR = x[num_comps*(i+1)];
    double uR   = x[num_comps*(i+1)+1] / rhoR;
    double ER   = x[num_comps*(i+1)+2] / rhoR;
    double pR = (gamma - 1.0) * rhoR * (ER - 0.5 * rhoR * uR * uR);

    // Solve the Riemann problem at the left and right faces to obtain
    // values for the fluxes.
    double rho, u, p;
    solve_riemann(gamma, rhoL, uL, pL, rhoR, uR, pR, &rho, &u, &p);

    // Now compute the inter-face flux.
    fluxes[3*i]   = rho*u;
    fluxes[3*i+1] = rho*u*u + p;
    fluxes[3*i+2] = rho*u*(0.5*u*u + p);
  }

  // Boundary fluxes are zero.
  fluxes[0] = fluxes[1] = fluxes[2] = 0.0;
  fluxes[3*N] = fluxes[3*N+1] = fluxes[3*N+2] = 0.0;

  // Compute the time derivatives of the conserved quantities using 
  // the discrete divergence theorem.
  for (int i = 0; i < N; ++i)
  {
    R[3*i]   = fluxes[3*i] - fluxes[3*(i+1)]/dx;
    R[3*i+1] = fluxes[3*i+1] - fluxes[3*(i+1)+1]/dx;
    R[3*i+2] = fluxes[3*i+2] - fluxes[3*(i+1)+2]/dx;
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
  return nonlinear_function_new("Sod shock tube", sod, 3, vtable); 
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

  // The adjacency graph for a 1D grid is fairly simple.
  adj_graph_t* graph = adj_graph_new(MPI_COMM_WORLD, N);
  for (int v = 0; v < adj_graph_num_vertices(graph); ++v)
  {
    if (v == 0)
    {
      adj_graph_set_num_edges(graph, v, 1);
      int* e = adj_graph_edges(graph, v);
      e[0] = 1;
    }
    else if (v == (N-1))
    {
      adj_graph_set_num_edges(graph, v, 1);
      int* e = adj_graph_edges(graph, v);
      e[0] = N-2;
    }
    else
    {
      adj_graph_set_num_edges(graph, v, 2);
      int* e = adj_graph_edges(graph, v);
      e[0] = v - 1;
      e[1] = v + 1;
    }
  }

  // Create a solution array and initialize it.
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
  const UnitTest tests[] = 
  {
    unit_test(test_sod_shock_tube)
  };
  return run_tests(tests);
}
