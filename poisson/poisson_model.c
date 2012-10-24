// poisson.c - The Poisson engine for Arbi.

#include <string.h>
#include <stdlib.h>
#include "poisson/poisson_model.h"
#include "core/constant_st_func.h"
#include "geometry/create_cvt.h"
#include "geometry/plane.h"
#include "geometry/cylinder.h"
#include "geometry/sphere.h"
#include "geometry/intersection.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct 
{
  mesh_t* mesh;         // Mesh.
  st_func_t* RHS;       // Right-hand side function.
} poisson_t;

// Benchmarks

static void run_paraboloid(options_t* opts)
{
  // Extract the dimension of the benchmark.
  char* dim_str = options_value(opts, "dim");
  int dim = 3;
  if (dim_str != NULL)
  {
    dim = atoi(options_value(opts, "dim"));
    if ((dim < 1) || (dim > 3))
      arbi_error("Invalid dimension: %d", dim);
  }

  // Get the geometry.
  static const char* geom_default = "dirichlet";
  char* geom = options_value(opts, "geometry");
  if (geom != NULL)
  {
    if (dim == 1)
    {
      arbi_warn("geometry is ignored for dim == 1");
    }
    else
    {
      if (strcmp(geom, "box") && 
          strcmp(geom, "cylinder") && 
          strcmp(geom, "sphere"))
      {
        arbi_error("Invalid geometry: %s", geom);
      }
    }
  }
  else
    geom = (char*)geom_default;

  // Get the boundary condition.
  static const char* bc_default = "dirichlet";
  char* bcond = options_value(opts, "bc");
  if (bcond != NULL)
  {
    if (strcmp(bcond, "dirichlet") && 
        strcmp(bcond, "neumann"))
    {
      arbi_error("Invalid boundary condition: %s", bcond);
    }
  }
  else
    bcond = (char*)bc_default;

  // Resolution.
  int N = 100;

  // Create the model.
  model_t* model = poisson_model_new(opts);
  poisson_t* p = model_context(model);

  // Create the (constant density) mesh.
  sp_func_t* boundary = NULL;
  if (!strcmp(geom, "box"))
  {
    sp_func_t* planes[6];
    vector_t n = {-1.0, 0.0, 0.0};
    point_t x = { 0.0, 0.5, 0.5};
    planes[0] = plane_new(n, x);
    n.x = 1.0, n.y = 0.0, n.z = 0.0;
    x.x = 1.0, x.y = 0.5, x.z = 0.5;
    planes[1] = plane_new(n, x);
    n.x = 0.0, n.y = -1.0, n.z = 0.0;
    x.x = 0.5, x.y = 0.0, x.z = 0.5;
    planes[2] = plane_new(n, x);
    n.x = 0.0, n.y = 1.0, n.z = 0.0;
    x.x = 0.5, x.y = 1.0, x.z = 0.5;
    planes[3] = plane_new(n, x);
    n.x = 0.0, n.y = 0.0, n.z = -1.0;
    x.x = 0.5, x.y = 0.5, x.z = 0.0;
    planes[4] = plane_new(n, x);
    n.x = 0.0, n.y = 0.0, n.z = 1.0;
    x.x = 0.5, x.y = 0.5, x.z = 1.0;
    planes[5] = plane_new(n, x);
    boundary = intersection_new(planes, 6);
  }
  else if (!strcmp(geom, "cylinder"))
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

  // Create the RHS function.
  double two = 2.0;
  p->RHS = constant_st_func_new(1, &two);

  // Set boundary conditions.
  // FIXME
  if (!strcmp(bcond, "dirichlet"))
  {
  }
  else if (!strcmp(bcond, "neumann"))
  {
  }

  // Run the thing.
  double t1 = 0.0, t2 = 1.0;
  model_run(model, t1, t2);

  // Clean up.
  model_free(model);
}

// Vtable stuff

static void poisson_run_benchmark(const char* benchmark, options_t* opts)
{
  if (!strcmp(benchmark, "paraboloid"))
  {
    run_paraboloid(opts);
  }
  else
  {
    char err[1024];
    snprintf(err, 1024, "poisson: unknown benchmark: '%s'", benchmark);
    arbi_error(err);
  }
}

static void poisson_init(void* p, double t)
{
}

static void poisson_advance(void* p, double t, double dt)
{
}

static void poisson_plot(void* p, plot_interface_t* plot, double t, int step)
{
}

static void poisson_dtor(void* p)
{
  free(p);
}

model_t* poisson_model_new(options_t* options)
{
  model_vtable vtable = { .run_benchmark = &poisson_run_benchmark,
                          .init = &poisson_init,
                          .advance = &poisson_advance,
                          .plot = &poisson_plot,
                          .dtor = &poisson_dtor};
  poisson_t* context = malloc(sizeof(poisson_t));
  return model_new("poisson", context, vtable);
}


#ifdef __cplusplus
}
#endif

