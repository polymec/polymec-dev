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
  model_compute_error_norms(model, solution, t2, lp_norms);
}

static void advect_run_1d_flow(options_t* options, 
                               st_func_t* velocity, 
                               st_func_t* diffusivity, 
                               st_func_t* source, 
                               st_func_t* initial_cond, 
                               st_func_t* solution, 
                               int dim)
{
  // Boundary conditions.
  str_ptr_unordered_map_t* bcs = str_ptr_unordered_map_new();

  // Dirichlet on -x/+x.
  str_ptr_unordered_map_insert(bcs, "-x", advect_bc_new(1.0, 0.0, solution));
  str_ptr_unordered_map_insert(bcs, "+x", advect_bc_new(1.0, 0.0, solution));

  // Transverse faces - homogeneous Neumann BCs.
  double z = 0.0;
  st_func_t* zero = constant_st_func_new(1, &z);
  str_ptr_unordered_map_insert(bcs, "-y", advect_bc_new(0.0, 1.0, zero));
  str_ptr_unordered_map_insert(bcs, "+y", advect_bc_new(0.0, 1.0, zero));
  str_ptr_unordered_map_insert(bcs, "-z", advect_bc_new(0.0, 1.0, zero));
  str_ptr_unordered_map_insert(bcs, "+z", advect_bc_new(0.0, 1.0, zero));

  // Run times.
  double t1 = 0.0, t2 = 1.0;

  // Base resolution, number of runs.
  int N0;
  int num_runs;
  switch(dim)
  {
    case 1: 
      N0 = 32;
      num_runs = 4;
      break;
    case 2:
      N0 = 16;
      num_runs = 3;
      break;
    case 3:
      N0 = 8;
      num_runs = 3;
      break;
  }

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
  int pos = 0;
  char* key;
  void* value;
  while (str_ptr_unordered_map_next(bcs, &pos, &key, &value))
    advect_bc_free(value);
  str_ptr_unordered_map_free(bcs);
  zero = NULL;
}

static void run_stationary_flow_1d(options_t* options)
{
  double z = 0.0, o = 1.0;
  double v[] = {1.0, 0.0, 0.0};
  st_func_t* zero = constant_st_func_new(1, &z);
  st_func_t* one = constant_st_func_new(1, &o);
  st_func_t* v0 = constant_st_func_new(3, v);
  advect_run_1d_flow(options, v0, zero, zero, one, one, 1);
  zero = one = v0 = NULL;
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
  advect_run_1d_flow(options, v0, one, zero, soln, soln, 1);
  zero = one = v0 = soln = NULL;
}

static void square_wave_1d_soln(void* ctx, point_t* x, double t, double* phi)
{
  double width = 0.25;
  double a = 1.0;
  *phi = (fabs(x->x - a*t) < 0.5*width) ? 1.0 : 0.0;
}

static void run_square_wave_1d(options_t* options)
{
  double z = 0.0;
  double v[] = {1.0, 0.0, 0.0};
  st_func_t* zero = constant_st_func_new(1, &z);
  st_func_t* v0 = constant_st_func_new(3, v);
  st_func_t* soln = st_func_from_func("square wave 1d", square_wave_1d_soln,
                                      ST_INHOMOGENEOUS, ST_NONCONSTANT, 1);
  advect_run_1d_flow(options, v0, zero, zero, soln, soln, 1);
  zero = v0 = soln = NULL;
}

void register_advect_benchmarks(model_t* model)
{
  model_register_benchmark(model, "stationary_flow_1d", run_stationary_flow_1d, "Stational flow in 1D (v = 1).");
  model_register_benchmark(model, "stationary_blayer_1d", run_stationary_blayer_1d, "Stational flow with boundary layer in 1D.");
  model_register_benchmark(model, "square_wave_1d", run_square_wave_1d, "Square wave propagation in 1D.");
}

#ifdef __cplusplus
}
#endif

