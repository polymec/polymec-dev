#include <string.h>
#include "core/unordered_map.h"
#include "core/constant_st_func.h"
#include "geometry/cubic_lattice.h"
#include "geometry/create_cubic_lattice_mesh.h"
#include "poisson/poisson_model.h"
#include "poisson/poisson_bc.h"
#include "poisson/interpreter_register_poisson_functions.h"
#include "poisson/register_poisson_benchmarks.h"

#ifdef __cplusplus
extern "C" {
#endif

static void run_analytic_problem(mesh_t* mesh, st_func_t* rhs, str_ptr_unordered_map_t* bcs, options_t* options, double t1, double t2, st_func_t* solution, double* lp_norms)
{
  // Create the model.
  model_t* model = create_poisson(mesh, rhs, bcs, solution, options);

  // Run the thing.
  model_run(model, t1, t2, INT_MAX);

  // Calculate the Lp norms of the error and write it to Lp_norms.
  model_compute_error_norms(model, solution, t2, lp_norms);

  // Clean up.
//  lp_norm = NULL;
  model_free(model);
}

static void laplace_1d_solution(void* context, point_t* x, double t, double* phi)
{
  phi[0] = 1.0 + 2.0*x->x;
}

static void laplace_1d_solution_grad(void* context, point_t* x, double t, double* grad_phi)
{
  grad_phi[0] = 2.0;
  grad_phi[1] = 0.0;
  grad_phi[2] = 0.0;
}

static void poisson_run_laplace_1d(options_t* options, int dim)
{
  // Parse any benchmark-specific options.
  bool all_dirichlet = false;
  bool reversed_bcs = false;
  char* bcs_opt = options_value(options, "bcs");
  if (bcs_opt != NULL)
  {
    if (!strcmp(bcs_opt, "dirichlet"))
      all_dirichlet = true;
    else if (!strcmp(bcs_opt, "reversed"))
      reversed_bcs = true;
  }

  // RHS function is zero for Laplace's equation.
  double z = 0.0;
  st_func_t* zero = constant_st_func_new(1, &z);

  // Analytic solution and gradient.
  st_func_t* sol = st_func_from_func("laplace_1d_sol", laplace_1d_solution,
                                     ST_INHOMOGENEOUS, ST_CONSTANT, 1);
  st_func_t* sol_grad = st_func_from_func("laplace_1d_sol_grad", laplace_1d_solution_grad,
                                          ST_INHOMOGENEOUS, ST_CONSTANT, 3);

  // Boundary conditions: Dirichlet on -x/+x (unless they've been reversed).
  str_ptr_unordered_map_t* bcs = str_ptr_unordered_map_new();
  if (!reversed_bcs)
  {
    str_ptr_unordered_map_insert(bcs, "-x", poisson_bc_new(1.0, 0.0, sol));
    str_ptr_unordered_map_insert(bcs, "+x", poisson_bc_new(1.0, 0.0, sol));
  }
  else
  {
    str_ptr_unordered_map_insert(bcs, "-x", poisson_bc_new(0.0, 1.0, sol_grad));
    str_ptr_unordered_map_insert(bcs, "+x", poisson_bc_new(0.0, 1.0, sol_grad));
  }

  // Transverse faces.
  if (all_dirichlet || reversed_bcs)
  {
    // Dirichlet BCs.
    str_ptr_unordered_map_insert(bcs, "-y", poisson_bc_new(1.0, 0.0, sol));
    str_ptr_unordered_map_insert(bcs, "+y", poisson_bc_new(1.0, 0.0, sol));
    str_ptr_unordered_map_insert(bcs, "-z", poisson_bc_new(1.0, 0.0, sol));
    str_ptr_unordered_map_insert(bcs, "+z", poisson_bc_new(1.0, 0.0, sol));
  }
  else
  {
    // Homogeneous Neumann BCs.
    str_ptr_unordered_map_insert(bcs, "-y", poisson_bc_new(0.0, 1.0, zero));
    str_ptr_unordered_map_insert(bcs, "+y", poisson_bc_new(0.0, 1.0, zero));
    str_ptr_unordered_map_insert(bcs, "-z", poisson_bc_new(0.0, 1.0, zero));
    str_ptr_unordered_map_insert(bcs, "+z", poisson_bc_new(0.0, 1.0, zero));
  }

  // Run time.
  double t = 0.0;

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
    run_analytic_problem(mesh, zero, bcs_copy, options, t, t, sol, Lp_norms[iter]);

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
    poisson_bc_free(value);
  str_ptr_unordered_map_free(bcs);
  zero = sol = sol_grad = NULL;
}

static void run_laplace_1d(options_t* options)
{
  poisson_run_laplace_1d(options, 1);
}

static void run_laplace_1d_2(options_t* options)
{
  poisson_run_laplace_1d(options, 2);
}

static void run_laplace_1d_3(options_t* options)
{
  poisson_run_laplace_1d(options, 3);
}

static void paraboloid_solution(void* context, point_t* x, double t, double* phi)
{
  double r2 = x->x*x->x + x->y*x->y; // Distance from center axis.
  phi[0] = 1.0 + r2;
//  printf("phi(%g, %g) = %g\n", x->x, x->y, phi[0]);
}

static void poisson_run_paraboloid(options_t* options, int dim)
{
  ASSERT((dim == 2) || (dim == 3));

  // Parse model-specific options.
  bool offcenter = false;
  bool all_dirichlet = false;

  char* bcs_opt = options_value(options, "bcs");
  if (bcs_opt != NULL)
  {
    if (!strcmp(bcs_opt, "dirichlet"))
      all_dirichlet = true;
  }
  char *geom = options_value(options, "geometry");
  if (geom != NULL)
  {
    // The mesh can be generated off-center so that the origin is
    // at the lower left.
    if (!strcasecmp(geom, "offcenter"))
      offcenter = true;
  }

  // RHS function.
  double four = 4.0;
  st_func_t* rhs = constant_st_func_new(1, &four);

  // Analytic solution.
  st_func_t* sol = st_func_from_func("paraboloid", paraboloid_solution,
                                     ST_INHOMOGENEOUS, ST_CONSTANT, 1);

  // Set up a Dirichlet boundary condition along each of the outside faces.
  str_ptr_unordered_map_t* bcs = str_ptr_unordered_map_new();
  str_ptr_unordered_map_insert(bcs, "+x", poisson_bc_new(1.0, 0.0, sol));
  str_ptr_unordered_map_insert(bcs, "-x", poisson_bc_new(1.0, 0.0, sol));
  str_ptr_unordered_map_insert(bcs, "+y", poisson_bc_new(1.0, 0.0, sol));
  str_ptr_unordered_map_insert(bcs, "-y", poisson_bc_new(1.0, 0.0, sol));
  
  double z = 0.0;
  st_func_t* zero = constant_st_func_new(1, &z);
  if (all_dirichlet)
  {
    str_ptr_unordered_map_insert(bcs, "+z", poisson_bc_new(1.0, 0.0, sol));
    str_ptr_unordered_map_insert(bcs, "-z", poisson_bc_new(1.0, 0.0, sol));
  }
  else
  {
    // Set up a homogeneous Neumann boundary condition on +/- z.
    str_ptr_unordered_map_insert(bcs, "+z", poisson_bc_new(0.0, 1.0, zero));
    str_ptr_unordered_map_insert(bcs, "-z", poisson_bc_new(0.0, 1.0, zero));
  }

  // Start/end time.
  double t = 0.0;

  // Base resolution and number of runs.
  int N0;
  int num_runs;
  switch (dim)
  {
    case 2:
      N0 = 16;
      num_runs = 4;
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
    if (dim > 1) Ny = Nx;
    if (dim > 2) Nz = Nx;
    bbox_t bbox;
    if (offcenter)
    {
      bbox.x1 = 0.0, bbox.x2 = 0.5, 
      bbox.y1 = 0.0, bbox.y2 = 0.5;
    }
    else
    {
      bbox.x1 = -0.5, bbox.x2 = 0.5, 
      bbox.y1 = -0.5, bbox.y2 = 0.5; 
    }
    bbox.z1 = 0.0, bbox.z2 = 1.0;
    if (dim == 2)
      bbox.z2 = 1.0/Nx;
    mesh_t* mesh = create_cubic_lattice_mesh_with_bbox(Nx, Ny, Nz, &bbox);
    tag_cubic_lattice_mesh_faces(mesh, Nx, Ny, Nz, "-x", "+x", "-y", "+y", "-z", "+z");
    str_ptr_unordered_map_t* bcs_copy = str_ptr_unordered_map_copy(bcs);
    run_analytic_problem(mesh, rhs, bcs_copy, options, t, t, sol, Lp_norms[iter]);

    // If we run in 2D, we need to adjust the norms.
    if (dim == 2)
    {
      Lp_norms[iter][1] *= Nx;
      Lp_norms[iter][2] *= Nx;
    }
    log_info("iteration %d (Nx = Ny = %d): L1 = %g, L2 = %g, Linf = %g", iter, Nx, Lp_norms[iter][1], Lp_norms[iter][2], Lp_norms[iter][0]);
  }

  // Clean up.
  int pos = 0;
  char* key;
  void* value;
  while (str_ptr_unordered_map_next(bcs, &pos, &key, &value))
    poisson_bc_free(value);
  str_ptr_unordered_map_free(bcs);
  zero = sol = rhs = NULL;
}

static void run_paraboloid(options_t* options)
{
  poisson_run_paraboloid(options, 2);
}

static void run_paraboloid_3(options_t* options)
{
  poisson_run_paraboloid(options, 3);
}

void register_poisson_benchmarks(model_t* model)
{
  model_register_benchmark(model, "laplace_1d", run_laplace_1d, "Laplace's equation in 1D Cartesian coordinates.");
  model_register_benchmark(model, "laplace_1d_2", run_laplace_1d_2, "Laplace's equation in 1D Cartesian coordinates (run in 2D).");
  model_register_benchmark(model, "laplace_1d_3", run_laplace_1d_3, "Laplace's equation in 1D Cartesian coordinates (run in 3D).");
  model_register_benchmark(model, "paraboloid", run_paraboloid, "A paraboloid solution to Poisson's equation (2D).");
  model_register_benchmark(model, "paraboloid_3", run_paraboloid_3, "A paraboloid solution to Poisson's equation (3D).");
}

#ifdef __cplusplus
}
#endif

