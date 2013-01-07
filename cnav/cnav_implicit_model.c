#include <string.h>
#include "core/unordered_map.h"
#include "io/silo_io.h"
#include "io/vtk_plot_io.h"
#include "io/gnuplot_io.h"
#include "cnav/cnav_implicit_model.h"
#include "cnav/cnav_bc.h"
#include "cnav/interpreter_register_cnav_functions.h"
#include "cnav/register_cnav_benchmarks.h"

#ifdef __cplusplus
extern "C" {
#endif

// This data structure represents all implicitly-integrated implementations
// of the compressible Navier-Stokes model.
typedef struct 
{
} cnav_implicit_t;

static double cnav_max_dt(void* context, double t, char* reason)
{
  cnav_implicit_t* cnav = (cnav_implicit_t*)context;
  return model_max_dt(cnav->model, reason);
}

static void cnav_advance(void* context, double t, double dt)
{
  cnav_implicit_t* cnav = (cnav_implicit_t*)context;
  model_advance(cnav->model, dt);
}

static void cnav_init(void* context, double t)
{
  cnav_implicit_t* cnav = (cnav_implicit_t*)context;
  model_init(cnav->model, t);
}

static void cnav_save(void* context, io_interface_t* io, double t, int step)
{
  cnav_implicit_t* cnav = (cnav_implicit_t*)context;
  model_save(cnav->model);
}

static void cnav_dtor(void* ctx)
{
  cnav_implicit_t* cnav = (cnav_implicit_t*)context;
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
                          .dtor = cnav_dtor};
  cnav_implicit_t* cnav = malloc(sizeof(cnav_implicit_t));
  model_t* model = model_new("Implicit compressible Navier-Stokes", cnav, vtable, options);

  // Set up the saver.
  io_interface_t* saver = silo_io_new(cnav->comm, 0, false);
  model_set_saver(model, saver);

  return model;
}

model_t* create_cnav_implicit(cnav_integrator_t integrator,
                              mesh_t* mesh,
                              cnav_eos_t* equation_of_state,
//                             reaction_network_t* reactions,
                              st_func_t* source, 
                              st_func_t* initial_cond, 
                              str_ptr_unordered_map_t* bcs, 
                              st_func_t* solution,
                              options_t* options)
{
  ASSERT(integrator != CNAV_SEMI_IMPLICIT);
  ASSERT(mesh != NULL);
  ASSERT(equation_of_state != NULL);
  ASSERT(source != NULL);
  ASSERT(initial_cond != NULL);
  ASSERT(st_func_num_comp(source) == st_func_num_comp(initial_cond));

  // Create the model.
  model_t* model = cnav_implicit_model_new(integrator, options);
  cnav_implicit_t* cnav = (cnav_implicit_t*)model_context(model);
  // FIXME
  return model;
}

#ifdef __cplusplus
}
#endif

