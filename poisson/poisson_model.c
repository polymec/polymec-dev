// poisson.c - The Poisson engine for Arbi.

#include <string.h>
#include <stdlib.h>
#include "petscmat.h"
#include "petscvec.h"
#include "core/unordered_map.h"
#include "core/constant_st_func.h"
#include "geometry/create_cubic_lattice_mesh.h"
#include "geometry/create_cvt.h"
#include "geometry/plane.h"
#include "geometry/cylinder.h"
#include "geometry/sphere.h"
#include "geometry/intersection.h"
#include "poisson/poisson_model.h"

#ifdef __cplusplus
extern "C" {
#endif

// Boundary condition structure.
typedef enum { DIRICHLET, NEUMANN, ROBIN } poisson_bc_type;
typedef struct
{
  poisson_bc_type type;
  st_func_t* func;
} poisson_bc_t;

static poisson_bc_t* create_dirichlet_bc(st_func_t* func)
{
  ASSERT(func != NULL);
  poisson_bc_t* bc = malloc(sizeof(poisson_bc_t));
  bc->type = DIRICHLET;
  bc->func = func;
  return bc;
}

static poisson_bc_t* create_neumann_bc(st_func_t* func)
{
  ASSERT(func != NULL);
  poisson_bc_t* bc = malloc(sizeof(poisson_bc_t));
  bc->type = NEUMANN;
  bc->func = func;
  return bc;
}

static poisson_bc_t* create_robin_bc(st_func_t* func)
{
  ASSERT(func != NULL);
  poisson_bc_t* bc = malloc(sizeof(poisson_bc_t));
  bc->type = ROBIN;
  bc->func = func;
  return bc;
}

static void free_bc(void* bc)
{
  poisson_bc_t* pbc = (poisson_bc_t*)bc;
  pbc->func = NULL;
}

// Context structure.
typedef struct 
{
  mesh_t* mesh;             // Mesh.
  Mat L;                    // Laplacian.
  st_func_t* rhs;           // Right-hand side function.
  double* phi;              // Solution.

  str_ptr_unordered_map_t* bcs;  // Boundary conditions.
} poisson_t;

// A proper constructor.
static model_t* create_poisson(mesh_t* mesh, st_func_t* rhs, str_ptr_unordered_map_t* bcs)
{
  model_t* poisson = poisson_model_new(NULL);
  poisson_t* p = (poisson_t*)model_context(poisson);
  p->mesh = mesh;
  p->rhs = rhs;
  p->bcs = bcs;
  return poisson;
}

//------------------------------------------------------------------------
//                            Benchmarks
//------------------------------------------------------------------------

// Mesh generation for benchmark problems.
static mesh_t* create_cube_mesh(int dim, int N)
{
  // Get the characteristic resolution.
  int N3[3] = {N, 1, 1};
  for (int d = 0; d < dim; ++d)
    N3[d] = N;

  return create_cubic_lattice_mesh(N3[0], N3[1], N3[2]);
}

static void run_analytic_problem(mesh_t* mesh, st_func_t* rhs, str_ptr_unordered_map_t* bcs, double t1, double t2, sp_func_t* solution, double* Lpnorms)
{
  // Create the model.
  model_t* model = create_poisson(mesh, rhs, bcs);

  // Run the thing.
  model_run(model, t1, t2);

  // Calculate the Lp norm of the error and write it to Lpnorms.

  // Clean up.
  model_free(model);
}

static void poisson_run_laplace_sov(int variant)
{
  // Dimension.
  int dim;
  // FIXME

  // RHS function.
  st_func_t* rhs;
  // FIXME

  // Boundary conditions.
  str_ptr_unordered_map_t* bcs = str_ptr_unordered_map_new();
  if ((variant % 3) == 0)
    str_ptr_unordered_map_insert(bcs, "boundary", create_dirichlet_bc(NULL));
  else if ((variant % 3) == 1)
    str_ptr_unordered_map_insert(bcs, "boundary", create_neumann_bc(NULL));
  else // if ((variant % 3) == 2)
    str_ptr_unordered_map_insert(bcs, "boundary", create_robin_bc(NULL));

  // Analytic solution.
  sp_func_t* sol;
  // FIXME = create_paraboloid_solution(variant);

  // Start/end times.
  double t1, t2;
  // FIXME

  // Base resolution.
  int N0;
  // FIXME

  // Number of refinements.
  int num_refinements;
  // FIXME

  // Do a convergence study.
  double Lpnorms[num_refinements][3];
  for (int iter = 0; iter < num_refinements; ++iter)
  {
    int N = pow(N0, iter+1);
    mesh_t* mesh = create_cube_mesh(dim, N);
    run_analytic_problem(mesh, rhs, bcs, t1, t2, sol, Lpnorms[iter]);
  }
}

static void poisson_run_paraboloid(int variant)
{
  // Dimension.
  int dim;
  // FIXME

  // RHS function.
  st_func_t* rhs;
  // FIXME = create_paraboloid_rhs(variant);

  // Boundary conditions.
  str_ptr_unordered_map_t* bcs = str_ptr_unordered_map_new();
  if ((variant % 3) == 0)
    str_ptr_unordered_map_insert(bcs, "boundary", create_dirichlet_bc(NULL));
  else if ((variant % 3) == 1)
    str_ptr_unordered_map_insert(bcs, "boundary", create_neumann_bc(NULL));
  else // if ((variant % 3) == 2)
    str_ptr_unordered_map_insert(bcs, "boundary", create_robin_bc(NULL));

  // Analytic solution.
  sp_func_t* sol;
  // FIXME = create_paraboloid_solution(variant);

  // Start/end times.
  double t1, t2;
  // FIXME

  // Base resolution.
  int N0;
  // FIXME

  // Number of refinements.
  int num_refinements;
  // FIXME

  // Do a convergence study.
  double Lpnorms[num_refinements][3];
  for (int iter = 0; iter < num_refinements; ++iter)
  {
    int N = pow(N0, iter+1);
    mesh_t* mesh = create_cube_mesh(dim, N);
    run_analytic_problem(mesh, rhs, bcs, t1, t2, sol, Lpnorms[iter]);
  }
}

//------------------------------------------------------------------------
//                        Model implementation
//------------------------------------------------------------------------

static void poisson_run_benchmark(const char* benchmark)
{
  char* variant_str = NULL;
  if ((variant_str = strstr(benchmark, "laplace_sov")) != NULL)
  {
    int variant = atoi(variant_str + strlen("laplace_sov"));
    poisson_run_laplace_sov(variant);
  }
  else if ((variant_str = strstr(benchmark, "paraboloid")) != NULL)
  {
    int variant = atoi(variant_str + strlen("paraboloid"));
    poisson_run_paraboloid(variant);
  }
  else
  {
    char err[1024];
    snprintf(err, 1024, "poisson: unknown benchmark: '%s'", benchmark);
    arbi_error(err);
  }
}

static void poisson_advance(void* p, double t, double dt)
{
}

static void poisson_init(void* p, double t)
{
  // We simply solve the problem for t = 0.
  poisson_advance(p, t, 0.0);
}

static void poisson_plot(void* p, plot_interface_t* plot, double t, int step)
{
}

static void poisson_dump(void* context, io_interface_t* io, double t, int step)
{
  ASSERT(context != NULL);
  poisson_t* p = (poisson_t*)context;

  io_dataset_t* dataset = io_dataset_new("default", 1, 0);
  io_dataset_write_mesh(dataset, p->mesh);
  io_dataset_write_field(dataset, "phi", p->phi, 1, MESH_CELL);
  io_append_dataset(io, dataset);
}

static void poisson_dtor(void* ctx)
{
  poisson_t* p = (poisson_t*)ctx;
  mesh_free(p->mesh);
  int pos = 0;
  char* key;
  void* val;
  while (str_ptr_unordered_map_next(p->bcs, &pos, &key, &val))
  {
    free(key);
    free_bc(val);
  }
  free(p);
}

model_t* poisson_model_new(options_t* options)
{
  model_vtable vtable = { .run_benchmark = &poisson_run_benchmark,
                          .init = &poisson_init,
                          .advance = &poisson_advance,
                          .plot = &poisson_plot,
                          .dump = &poisson_dump,
                          .dtor = &poisson_dtor};
  poisson_t* context = malloc(sizeof(poisson_t));
  model_t* model = model_new("poisson", context, vtable);
  static const char* benchmarks[] = {"laplace_sov1", "laplace_sov2", "laplace_sov3", 
                                     "paraboloid1", "paraboloid2", "paraboloid3",
                                     NULL};
  model_register_benchmarks(model, benchmarks);
  if (options != NULL)
  {
  }
  return model;
}


#ifdef __cplusplus
}
#endif

