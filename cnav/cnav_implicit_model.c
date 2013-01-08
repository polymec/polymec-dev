#include <string.h>
#include "core/unordered_map.h"
#include "io/silo_io.h"
#include "io/vtk_plot_io.h"
#include "io/gnuplot_io.h"
#include "cnav/cnav_implicit_model.h"
#include "cnav/cnav_bc.h"
#include "petscts.h"

#ifdef __cplusplus
extern "C" {
#endif

// This data structure represents all implicitly-integrated implementations
// of the compressible Navier-Stokes model.
typedef struct 
{
  MPI_Comm comm;
  mesh_t* mesh;
  TS stepper;            // PETSC time stepper.
  Vec U;                 // Computed solution.
  str_ptr_unordered_map_t* bcs; // Boundary conditions.
  boundary_cell_map_t* boundary_cells; // Boundary cell data.
  bool initialized;
  double abs_tol, rel_tol;  // Absolute and relative tolerances.
  double CFL;               // CFL safety factor.
} cnav_implicit_t;

static PetscErrorCode compute_F(TS stepper, PetscReal t, Vec U, Vec U_dot, Vec F, void* context)
{
  // FIXME
  return 0;
}

static void cnav_advance(void* context, double t, double dt)
{
  cnav_implicit_t* cnav = (cnav_implicit_t*)context;
  // Take a single step.
  TSStep(cnav->stepper);
}

static void cnav_reconnect(void* context, mesh_t* new_mesh)
{
}

static void cnav_init(void* context, double t)
{
  cnav_implicit_t* cnav = (cnav_implicit_t*)context;

  // Make sure our solution vector is allocated.
  if (cnav->initialized)
  {
    VecDestroy(&cnav->U);
    cnav->initialized = false;
  }
  VecCreate(cnav->comm, &cnav->U);
  VecSetType(cnav->U, VECSEQ);
  VecSetSizes(cnav->U, cnav->mesh->num_cells, PETSC_DECIDE);

  // Initialize the solution.
  // FIXME
  double* U;
  VecGetArray(cnav->U, &U);
  for (int c = 0; c < num_cells; ++c)
    st_func_eval(a->initial_cond, &a->mesh->cells[c].center, t, &U[c]);
  VecRestoreArray(cnav->U, &U);
  TSSetSolution(cnav->stepper, cnav->U);

  // Find the initial time step size.
  double dt = 1e-3; // FIXME
  TSSetInitialTimeStep(cnav->stepper, t, dt); 

  cnav->initialized = true;
}

static void cnav_save(void* context, io_interface_t* io, double t, int step)
{
  cnav_implicit_t* cnav = (cnav_implicit_t*)context;
  model_save(cnav->model);
}

static void cnav_dtor(void* ctx)
{
  cnav_implicit_t* cnav = (cnav_implicit_t*)context;
  TSDestroy(cnav->stepper);
  VecDestroy(cnav->U);
  free(cnav);
}

model_t* cnav_implicit_model_new(cnav_integrator_t integrator,
                                 options_t* options)
{
  ASSERT(integrator != CNAV_SEMI_IMPLICIT);

  model_vtable vtable = { .read_input = cnav_read_input,
                          .init = cnav_init,
                          .advance = cnav_advance,
                          .save = cnav_save,
                          .dtor = cnav_dtor};
  cnav_implicit_t* cnav = malloc(sizeof(cnav_implicit_t));
  model_t* model = model_new("Implicit compressible Navier-Stokes", cnav, vtable, options);

  // Initialize the bookkeeping structures.
  cnav->comm = MPI_COMM_WORLD;
  TSCreate(cnav->comm, &cnav->stepper);
  TSSetProblemType(cnav->stepper, TS_NONLINEAR);
  TSSetIFFunction(cnav->stepper, compute_F, cnav);
  TSType int_type;
  if (integrator == IMPLICIT_BACKWARD_EULER)
    int_type = TSBEULER;
  else 
    int_type = TSSUNDIALS;
  TSSetType(cnav->stepper, int_type);
  if (int_type == TSSUNDIALS)
  {
    if (integrator == IMPLICIT_ADAMS)
      TSSundialsSetType(cnav->stepper, SUNDIALS_ADAMS);
    else
      TSSundialsSetType(cnav->stepper, SUNDIALS_BDF);
  }
    
  cnav->initialized = false;

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
  ASSERT(mesh != NULL);
  ASSERT(equation_of_state != NULL);
  ASSERT(source != NULL);
  ASSERT(initial_cond != NULL);
  ASSERT(st_func_num_comp(source) == st_func_num_comp(initial_cond));

  // Create the model.
  model_t* model = cnav_implicit_model_new(integrator, options);
  cnav_implicit_t* cnav = (cnav_implicit_t*)model_context(model);
  cnav->mesh = mesh;
  // FIXME
  return model;
}

#ifdef __cplusplus
}
#endif

