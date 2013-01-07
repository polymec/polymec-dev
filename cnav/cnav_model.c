#include <string.h>
#include "core/unordered_map.h"
#include "io/silo_io.h"
#include "io/vtk_plot_io.h"
#include "io/gnuplot_io.h"
#include "cnav/cnav_model.h"
#include "cnav/cnav_bc.h"
#include "cnav/interpreter_register_cnav_functions.h"
#include "cnav/register_cnav_benchmarks.h"

#ifdef __cplusplus
extern "C" {
#endif

// The generic "cnav" model is a wrapper around different implementations
// of the compressible Navier-Stokes model.
typedef struct 
{
  model_t* model;                    // The actual model.
} cnav_t;

static double cnav_max_dt(void* context, double t, char* reason)
{
  cnav_t* cnav = (cnav_t*)context;
  return model_max_dt(cnav->model, reason);
}

static void cnav_advance(void* context, double t, double dt)
{
  cnav_t* cnav = (cnav_t*)context;
  model_advance(cnav->model, dt);
}

static void cnav_read_input(void* context, interpreter_t* interp, options_t* options)
{
  cnav_t* cnav = (cnav_t*)context;
  // FIXME: Figure out how to pass this through!
}

static void cnav_init(void* context, double t)
{
  cnav_t* cnav = (cnav_t*)context;
  model_init(cnav->model, t);
}

static void cnav_plot(void* context, io_interface_t* io, double t, int step)
{
  cnav_t* cnav = (cnav_t*)context;
  // FIXME: Actual stuff goes here!
}

static void cnav_save(void* context, io_interface_t* io, double t, int step)
{
  cnav_t* cnav = (cnav_t*)context;
  model_save(cnav->model);
}

static void cnav_compute_error_norms(void* context, st_func_t* solution, double t, double* lp_norms)
{
  cnav_t* cnav = (cnav_t*)context;
  // FIXME: Actual stuff goes here!
}

static void cnav_dtor(void* ctx)
{
  cnav_t* cnav = (cnav_t*)context;
  model_free(cnav->model);
  free(cnav);
}

model_t* cnav_model_new(options_t* options)
{
  model_vtable vtable = { .read_input = cnav_read_input,
                          .init = cnav_init,
                          .max_dt = cnav_max_dt,
                          .advance = cnav_advance,
                          .save = cnav_save,
                          .plot = cnav_plot,
                          .compute_error_norms = cnav_compute_error_norms,
                          .dtor = cnav_dtor};
  cnav_t* cnav = malloc(sizeof(cnav_t));
  cnav->model = NULL;

  // Register benchmarks.
  register_cnav_benchmarks(model);

  // Set up the plotter.
  io_interface_t* plotter = NULL;
  char* which_plotter = options_value(options, "plotter");
  if (which_plotter != NULL)
  {
    if (!strcasecmp(which_plotter, "vtk"))
      plotter = vtk_plot_io_new(a->comm, 0, false);
    else if (!strcasecmp(which_plotter, "silo"))
      plotter = silo_plot_io_new(a->comm, 0, false);
    else if (!strcasecmp(which_plotter, "gnuplot"))
      plotter = gnuplot_io_new();
  }
  else
    plotter = vtk_plot_io_new(a->comm, 0, false);
  if (plotter != NULL)
  {
    log_detail("Setting plotter to '%s'...", which_plotter);
    model_set_plotter(model, plotter);
  }

  return model;
}

model_t* create_cnav(cnav_integrator_t integrator,
                     mesh_t* mesh,
                     cnav_eos_t* equation_of_state,
//                     reaction_network_t* reactions,
                     st_func_t* source, 
                     st_func_t* initial_cond, 
                     str_ptr_unordered_map_t* bcs, 
                     st_func_t* solution,
                     options_t* options)
{
  ASSERT(mesh != NULL);
  ASSERT(equation_of_state != NULL);
  ASSERT(source != NULL);
  ASSERT(initial_cond != NULL);
  ASSERT(st_func_num_comp(source) == st_func_num_comp(initial_cond));

  // Create the model.
  model_t* model = cnav_model_new(options);
  cnav_t* cnav = (cnav_t*)model_context(model);
  if (integrator == CNAV_SEMI_IMPLICIT)
    cnav->model = create_cnav_semi_implicit(mesh, equation_of_state, source, initial_cond, bcs, solution, options);
  else
    cnav->model = create_cnav_implicit(integrator, mesh, equation_of_state, source, initial_cond, bcs, solution, options);
  return model;
}

#ifdef __cplusplus
}
#endif

