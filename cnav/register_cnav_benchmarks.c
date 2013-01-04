#include <string.h>
#include "core/unordered_map.h"
#include "core/constant_st_func.h"
#include "geometry/cubic_lattice.h"
#include "geometry/create_cubic_lattice_mesh.h"
#include "advect/advect_model.h"
#include "advect/advect_bc.h"

#ifdef __cplusplus
extern "C" {
#endif

static void run_analytic_problem(model_t* model, 
                                 double t1, 
                                 double t2, 
                                 st_func_t* solution, 
                                 options_t* options,
                                 double* lp_norms)
{
  int max_steps = INT_MAX;
  char* max_steps_s = options_value(options, "max_steps");
  if (max_steps_s != NULL)
    max_steps = atoi(max_steps_s);

  // Run the thing.
  model_run(model, t1, t2, max_steps);

  // Calculate the Lp norms of the error and write it to Lp_norms.
  model_compute_error_norms(model, solution, lp_norms);
}

static void advect_run_1d_flow(options_t* options, 
                               st_func_t* velocity, 
                               st_func_t* diffusivity, 
                               st_func_t* source, 
                               st_func_t* initial_cond, 
                               void* left_bc,
                               void* right_bc,
                               st_func_t* solution, 
                               double t1,
                               double t2,
                               int dim, 
                               int num_runs)
{
  // Boundary conditions.
  str_ptr_unordered_map_t* bcs = str_ptr_unordered_map_new();

  // Dirichlet on -x/+x.
  str_ptr_unordered_map_insert(bcs, "-x", left_bc);
  str_ptr_unordered_map_insert(bcs, "+x", right_bc);

  // Transverse faces - homogeneous Neumann BCs.
  double z = 0.0;
  st_func_t* zero = constant_st_func_new(1, &z);
  str_ptr_unordered_map_insert(bcs, "-y", advect_bc_new(0.0, 1.0, zero));
  str_ptr_unordered_map_insert(bcs, "+y", advect_bc_new(0.0, 1.0, zero));
  str_ptr_unordered_map_insert(bcs, "-z", advect_bc_new(0.0, 1.0, zero));
  str_ptr_unordered_map_insert(bcs, "+z", advect_bc_new(0.0, 1.0, zero));

  // Base resolution.
  int N0;
  switch(dim)
  {
    case 1: 
      N0 = 32;
      break;
    case 2:
      N0 = 16;
      break;
    case 3:
      N0 = 8;
      break;
  }

  // Override num_runs?
  char* num_runs_str = options_value(options, "num_runs");
  if (num_runs_str != NULL)
    num_runs = atoi(num_runs_str);

  // Do a convergence study.
  double Lp_norms[num_runs][3];
  for (int iter = 0; iter < num_runs; ++iter)
  {
    int Nx = N0 * pow(2, iter), Ny = 1, Nz = 1;
    double dx = 1.0/Nx;
    bbox_t bbox = {.x1 = 0.0, .x2 = 1.0, .y1 = 0.0, .y2 = 1.0, .z1 = 0.0, .z2 = 1.0};
    if (dim == 1)
      bbox.y2 = bbox.z2 = dx;
    if (dim == 2)
    {
      Ny = Nx;
      bbox.z2 = dx;
    }
    if (dim == 3)
      Nz = Nx;
    mesh_t* mesh = create_cubic_lattice_mesh_with_bbox(Nx, Ny, Nz, &bbox);
    tag_cubic_lattice_mesh_faces(mesh, Nx, Ny, Nz, "-x", "+x", "-y", "+y", "-z", "+z");
    str_ptr_unordered_map_t* bcs_copy = str_ptr_unordered_map_copy(bcs);

    model_t* model = create_advect(mesh, velocity, diffusivity, source, 
                                   initial_cond, bcs_copy, solution, options);
    run_analytic_problem(model, t1, t2, solution, options, Lp_norms[iter]);
    model_free(model);

    // If we run in 1D or 2D, we need to adjust the norms.
    if (dim == 1)
    {
      Lp_norms[iter][1] *= Nx*Nx;
      Lp_norms[iter][2] *= Nx*Nx;
    }
    else if (dim == 2)
    {
      Lp_norms[iter][1] *= Nx;
      Lp_norms[iter][2] *= Nx;
    }
    log_urgent("iteration %d (Nx = %d): L1 = %g, L2 = %g, Linf = %g", iter, Nx, Lp_norms[iter][1], Lp_norms[iter][2], Lp_norms[iter][0]);
  }

  // Clean up.
  str_ptr_unordered_map_free(bcs);
}

static void run_stationary_flow_1d(options_t* options)
{
  double z = 0.0, o = 1.0;
  double v[] = {1.0, 0.0, 0.0};
  st_func_t* zero = constant_st_func_new(1, &z);
  st_func_t* one = constant_st_func_new(1, &o);
  st_func_t* v0 = constant_st_func_new(3, v);
  advect_bc_t* bc = advect_bc_new(1.0, 0.0, one);
  advect_run_1d_flow(options, v0, zero, zero, one, bc, bc, one, 0.0, 1.0, 1, 4);
}

static void stationary_blayer_1d_soln(void* ctx, point_t* x, double t, double* phi)
{
  double xx = x->x;
  double a = 1.0, d = 1.0;
  *phi = (exp(a/d) - exp(a*xx/d)) / (exp(a/d) - 1.0);
}

static void run_stationary_blayer_1d(options_t* options)
{
  double z = 0.0, o = 1.0;
  double v[] = {1.0, 0.0, 0.0};
  st_func_t* zero = constant_st_func_new(1, &z);
  st_func_t* one = constant_st_func_new(1, &o);
  st_func_t* v0 = constant_st_func_new(3, v);
  st_func_t* soln = st_func_from_func("blayer 1d", stationary_blayer_1d_soln,
                                      ST_INHOMOGENEOUS, ST_CONSTANT, 1);
  advect_bc_t* bc = advect_bc_new(1.0, 0.0, one);
  advect_run_1d_flow(options, v0, one, zero, soln, bc, bc, soln, 0.0, 1.0, 1, 4);
}

static void square_wave_1d_soln(void* ctx, point_t* x, double t, double* phi)
{
  // Taken from Toro (2009).
  double width = 0.4;
  double a = 1.0;
  double y = x->x - a*t - 0.5; // Center of the wave.
  while (y < 0.0) y += 1.0; // account for periodicity.
  *phi = ((fabs(y) < 0.5*width) || (fabs(y - 1.0) < 0.5*width)) ? 1.0 : 0.0;
}

static void run_square_wave_1d(options_t* options)
{
  double z = 0.0;
  double v[] = {1.0, 0.0, 0.0};
  st_func_t* zero = constant_st_func_new(1, &z);
  st_func_t* v0 = constant_st_func_new(3, v);
  st_func_t* soln = st_func_from_func("square wave 1d", square_wave_1d_soln,
                                      ST_INHOMOGENEOUS, ST_NONCONSTANT, 1);
  periodic_bc_t* bc = cubic_lattice_x_periodic_bc_new("-x", "+x");
  advect_run_1d_flow(options, v0, zero, zero, soln, bc, bc, soln, 0.0, 10.0, 1, 4);
}

static void gaussian_wave_1d_soln(void* ctx, point_t* x, double t, double* phi)
{
  // Taken from Toro (2009).
  double alpha = 1.0;
  double beta = 16.0;
  double a = 1.0;
  double y = x->x - a*t - 0.5; // Center of the wave.
  while (y < -0.5) y += 1.0; // account for periodicity.
  *phi = alpha * exp(-beta*y*y);
}

static void run_gaussian_wave_1d(options_t* options)
{
  double z = 0.0;
  double v[] = {1.0, 0.0, 0.0};
  st_func_t* zero = constant_st_func_new(1, &z);
  st_func_t* v0 = constant_st_func_new(3, v);
  st_func_t* soln = st_func_from_func("gaussian wave 1d", gaussian_wave_1d_soln,
                                      ST_INHOMOGENEOUS, ST_NONCONSTANT, 1);
  periodic_bc_t* bc = cubic_lattice_x_periodic_bc_new("-x", "+x");
  advect_run_1d_flow(options, v0, zero, zero, soln, bc, bc, soln, 0.0, 10.0, 1, 4);
}

static void sine_wave_1d_soln(void* ctx, point_t* x, double t, double* phi)
{
  double a = 1.0;
  *phi = sin(2.0*M_PI*(x->x - a*t));
}

static void run_sine_wave_1d(options_t* options)
{
  double z = 0.0;
  double v[] = {1.0, 0.0, 0.0};
  st_func_t* zero = constant_st_func_new(1, &z);
  st_func_t* v0 = constant_st_func_new(3, v);
  st_func_t* soln = st_func_from_func("sine wave 1d", sine_wave_1d_soln,
                                      ST_INHOMOGENEOUS, ST_NONCONSTANT, 1);
  periodic_bc_t* bc = cubic_lattice_x_periodic_bc_new("-x", "+x");
  advect_run_1d_flow(options, v0, zero, zero, soln, bc, bc, soln, 0.0, 5.0, 1, 4);
}

static void square_diffusion_1d_soln(void* ctx, point_t* x, double t, double* phi)
{
  double L = 1.0;  // Length of domain.
  double D = 0.01; // Diffusivity.
  double Vx = 0.0; // x velocity.
  double y = x->x - Vx*t - 0.5; // Center of the wave.
  double half_width = 0.25;     // half-width of the wave.
  if (t == 0.0)
  {
    *phi = (fabs(y) < half_width) ? 1.0 : 0.0;
  }
  else
  {
    double wL = half_width * L;
    double nu = D/(wL*wL);
    *phi = 0.5 * (erf((1.0 + fabs(y)/wL) / (2.0*sqrt(nu*t))) +  // left wave
                  erf((1.0 - fabs(y)/wL) / (2.0*sqrt(nu*t))));   // right wave
  }
}

static void run_square_diffusion_1d(options_t* options)
{
  double z = 0.0;
  double Dval = 0.01;
  double v[] = {0.0, 0.0, 0.0};

  char* D_str = options_value(options, "diffusivity");
  if (D_str != NULL)
    Dval = atof(D_str);

  st_func_t* zero = constant_st_func_new(1, &z);
  st_func_t* D = constant_st_func_new(1, &Dval);
  st_func_t* v0 = constant_st_func_new(3, v);
  st_func_t* soln = st_func_from_func("square diffusion 1d", square_diffusion_1d_soln,
                                      ST_INHOMOGENEOUS, ST_NONCONSTANT, 1);
  advect_bc_t* bc = advect_bc_new(1.0, 0.0, soln); // Dirichlet BC
  advect_run_1d_flow(options, v0, D, zero, soln, bc, bc, soln, 0.0, 2.0, 1, 4);
}

static void gaussian_diffusion_1d_soln(void* ctx, point_t* x, double t, double* phi)
{
  double D = 0.01; // Diffusivity.
  double t0 = 0.5; // Initial time.
  double Vx = 0.0; // x velocity.
  double y = x->x - Vx*t - 0.5; // Center of the wave.
  *phi = sqrt(t0/t) * exp(-y*y/(4.0*D*t));
}

static void run_gaussian_diffusion_1d(options_t* options)
{
  double z = 0.0;
  double Dval = 0.01;
  double v[] = {0.0, 0.0, 0.0};

  char* D_str = options_value(options, "diffusivity");
  if (D_str != NULL)
    Dval = atof(D_str);

  st_func_t* zero = constant_st_func_new(1, &z);
  st_func_t* D = constant_st_func_new(1, &Dval);
  st_func_t* v0 = constant_st_func_new(3, v);
  st_func_t* soln = st_func_from_func("gaussian diffusion 1d", gaussian_diffusion_1d_soln,
                                      ST_INHOMOGENEOUS, ST_NONCONSTANT, 1);
  advect_bc_t* bc = advect_bc_new(1.0, 0.0, soln); // Dirichlet BC
  advect_run_1d_flow(options, v0, D, zero, soln, bc, bc, soln, 0.5, 1.0, 1, 4);
}

void register_advect_benchmarks(model_t* model)
{
  model_register_benchmark(model, "stationary_flow_1d", run_stationary_flow_1d, "Stational flow in 1D (v = 1).");
  model_register_benchmark(model, "stationary_blayer_1d", run_stationary_blayer_1d, "Stational flow with boundary layer in 1D.");
  model_register_benchmark(model, "square_wave_1d", run_square_wave_1d, "Square wave propagation in 1D.");
  model_register_benchmark(model, "gaussian_wave_1d", run_gaussian_wave_1d, "Gaussian wave propagation in 1D.");
  model_register_benchmark(model, "sine_wave_1d", run_sine_wave_1d, "Sine wave propagation in 1D.");
  model_register_benchmark(model, "square_diffusion_1d", run_square_diffusion_1d, "Square wave diffusing in 1D.");
  model_register_benchmark(model, "gaussian_diffusion_1d", run_gaussian_diffusion_1d, "Gaussian wave diffusing in 1D.");
}

#ifdef __cplusplus
}
#endif

