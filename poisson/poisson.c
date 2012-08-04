// poisson.c - The Poisson engine for Arbi.

#include <stdlib.h>
#include "poisson/poisson.h"
#include "core/st_func.h"
//#include "core/lin_solver.h"

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
  int dim = atoi(options_value(opts, "dim"));
  if ((dim < 1) || (dim > 3))
    arbi_error("Invalid dimension: %d", dim);

  // Create the model.
  model_t* poisson = model_new("poisson", opts);

  // Clean up.
  model_free(poisson);
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

