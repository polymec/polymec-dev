// poisson.c - The Poisson engine for Arbi.

#include <string.h>
#include <stdlib.h>
#include "petscmat.h"
#include "petscvec.h"
#include "core/hash_map.h"
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

  str_ptr_hash_map_t* bcs;  // Boundary conditions.
} poisson_t;

// A proper constructor.
static model_t* create_poisson(mesh_t* mesh, st_func_t* rhs, str_ptr_hash_map_t* bcs)
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

// Paraboloid and variants.
static mesh_t* create_paraboloid_mesh(int variant, int iteration)
{
  int dim = (variant % 3) + 1;

  // Get the characteristic resolution.
  int N0 = (dim == 1) ? 128 : (dim == 2) ? 64 : 16;
  int N[3] = {N0, N0, N0};
  for (int iter = 0; iter < iteration; ++iter)
  {
    for (int d = 0; d < dim; ++d)
      N[d] *= 2;
  }

  bool is_box = (variant < 6);
  if (is_box)
  {
    return create_cubic_lattice_mesh(N[0], N[1], N[2]);
  }
  else
  {
#if 0
    sp_func_t* boundary = NULL;
    if (!strcmp(geom, "cylinder"))
    {
      sp_func_t* cyl[3];
      vector_t n = { 0.0, 0.0,-1.0};
      point_t x = { 0.5, 0.5, 0.0};
      cyl[0] = plane_new(n, x);
      n.x = 0.0, n.y = 0.0, n.z = 1.0;
      x.x = 0.5, x.y = 0.5, x.z = 1.0;
      cyl[1] = plane_new(n, x);
      cyl[2] = cylinder_new(n, x, 1.0);
      boundary = intersection_new(cyl, 3);
    }
    else if (!strcmp(geom, "sphere"))
    {
      point_t x = { 0.5, 0.5, 0.5};
      boundary = sphere_new(x, 1.0);
    }
    double one = 1.0;
    sp_func_t* rho = constant_sp_func_new(1, &one);
    p->mesh = create_bounded_cvt(N, rho, &cvt_energy_function, boundary);
#endif
    return NULL;
  }
}

static st_func_t* create_paraboloid_rhs(int variant)
{
  double one = 1.0;
  return constant_st_func_new(1, &one);
}

static str_ptr_hash_map_t* create_paraboloid_bcs(int variant)
{
  str_ptr_hash_map_t* bcs = str_ptr_hash_map_new();
  if ((variant % 3) == 0)
    str_ptr_hash_map_insert(bcs, "boundary", create_dirichlet_bc(NULL));
  else if ((variant % 3) == 1)
    str_ptr_hash_map_insert(bcs, "boundary", create_neumann_bc(NULL));
  else // if ((variant % 3) == 2)
    str_ptr_hash_map_insert(bcs, "boundary", create_robin_bc(NULL));
  return bcs;
}

static void set_paraboloid_times(int variant, double* t1, double* t2)
{
  *t1 = 0.0;
  *t2 = 1.0;
}

static sp_func_t* create_paraboloid_solution(int variant)
{
  return NULL;
}

static void run_analytic_problem(mesh_t* mesh, st_func_t* rhs, str_ptr_hash_map_t* bcs, double t1, double t2, sp_func_t* solution, double* Lpnorms)
{
  // Create the model.
  model_t* model = create_poisson(mesh, rhs, bcs);

  // Run the thing.
  model_run(model, t1, t2);

  // Calculate the Lp norm of the error and write it to Lpnorms.

  // Clean up.
  model_free(model);
}

//------------------------------------------------------------------------
//                        Model implementation
//------------------------------------------------------------------------

static void poisson_run_benchmark(const char* benchmark)
{
  char* variant_str = NULL;
  if ((variant_str = strstr(benchmark, "paraboloid")) != NULL)
  {
    double Lpnorms[4][3];
    int variant = atoi(variant_str + strlen("paraboloid"));
    st_func_t* rhs = create_paraboloid_rhs(variant);
    str_ptr_hash_map_t* bcs = create_paraboloid_bcs(variant);
    sp_func_t* sol = create_paraboloid_solution(variant);
    double t1, t2;
    set_paraboloid_times(variant, &t1, &t2);
    for (int iter = 0; iter < 4; ++iter)
    {
      mesh_t* mesh = create_paraboloid_mesh(variant, iter);
      run_analytic_problem(mesh, rhs, bcs, t1, t2, sol, Lpnorms[iter]);
    }
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
  // FIXME: BCs are leaked.
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
  if (options != NULL)
  {
  }
  return model;
}


#ifdef __cplusplus
}
#endif

