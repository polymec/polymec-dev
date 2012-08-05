// poisson.c - The Poisson engine for Arbi.

#include <string.h>
#include <stdlib.h>
#include "poisson/poisson.h"
#include "core/constant_st_func.h"
#include "geometry/create_box_mesh.h"
#include "geometry/create_cyl_mesh.h"
#include "geometry/create_sphere_mesh.h"

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

  // Create the model.
  model_t* model = model_new("poisson", opts);
  poisson_t* p = model_context(model);

  // Create the mesh.
  if ((dim == 1) || !strcmp(geom, "box"))
  {
    int N[3] = {10, 1, 1};
    if (dim >= 2)
      N[1] = 10;
    if (dim == 3)
      N[2] = 10;
    double low[3] = {0.0, 0.0, 0.0};
    double high[3] = {1.0, 1.0, 1.0};
    p->mesh = create_box_mesh(N, low, high);
  }
  else if (!strcmp(geom, "cylinder"))
  {
    int Ncenter = 10, Nradial = 10, Naxial = 10;
    if (dim == 2) 
      Naxial = 1;
    double Lbox = 0.5, R = 1.0, Lz = 1.0;
    p->mesh = create_cyl_mesh(Ncenter, Nradial, Naxial, Lbox, R, Lz);
  }
  else if (!strcmp(geom, "sphere"))
  {
    int Ncenter = 10, Nradial = 10;
    double Lbox = 0.5, R = 1.0;
    p->mesh = create_sphere_mesh(Ncenter, Nradial, Lbox, R);
  }

  // Create the RHS function.
  double rho = 2.0;
  p->RHS = constant_st_func_new(1, &rho);

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
static void* poisson_ctor(options_t* opts)
{
  poisson_t* p = malloc(sizeof(poisson_t));
  return p;
}

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

static model_vtable poisson_vtable = 
{
  .ctor          = &poisson_ctor,
  .run_benchmark = &poisson_run_benchmark,
  .init          = &poisson_init,
  .advance       = &poisson_advance,
  .plot          = &poisson_plot,
  .dtor          = &poisson_dtor
};

// FIXME
static const char* poisson_usage = "poisson model\n\n";

void register_poisson()
{
  register_model("poisson", poisson_usage, poisson_vtable);
}

#ifdef __cplusplus
}
#endif

