#include <string.h>

#include "cvode/cvode.h"
#if USE_MPI
#include "nvector/nvector_parallel.h"
#else
#include "nvector/nvector_serial.h"
#endif

#include "core/unordered_map.h"
#include "io/silo_io.h"
#include "io/vtk_plot_io.h"
#include "io/gnuplot_io.h"
#include "cnav/cnav_implicit_model.h"
#include "cnav/cnav_bc.h"

#ifdef __cplusplus
extern "C" {
#endif

// This data structure represents all implicitly-integrated implementations
// of the compressible Navier-Stokes model.
typedef struct 
{
#if USE_MPI
  MPI_Comm comm;
#endif
  mesh_t* mesh;
  void* cvode;                // Workspace for CVODE.
  N_Vector U;                 // Computed solution.
  str_ptr_unordered_map_t* bcs; // Boundary conditions.
  boundary_cell_map_t* boundary_cells; // Boundary cell data.
  double abs_tol, rel_tol;  // Absolute and relative tolerances.
  double CFL;               // CFL safety factor.
} cnav_implicit_t;

static inline N_Vector N_Vector_new(MPI_Comm comm, int dim)
{
#ifdef USE_MPI
  int tot; // FIXME
  return N_VNew_Parallel(comm, dim, tot);
#else
  return N_VNew_Serial(dim);
#endif
}

static inline N_Vector N_Vector_free(N_Vector* v)
{
#ifdef USE_MPI
  return N_VDestroy_Parallel(v);
#else
  return N_VDestroy_Serial(v);
#endif
}

static inline double* N_Vector_data(N_Vector v)
{
#ifdef USE_MPI
  return NV_DATA_P(v);
#else
  return NV_DATA_S(v);
#endif
}

static int compute_F(double t, N_Vector U, N_Vector U_dot, void* context)
{
  // FIXME
  return 0;
}

static double cnav_max_dt(void* context, double t, char* reason)
{
  cnav_implicit_t* cnav = (cnav_implicit_t*)context;

  // Find the minimum cell length / velocity ratio.
  double dt = FLT_MAX;
  // FIXME
  return dt;
}

static void cnav_advance(void* context, double t, double dt)
{
  cnav_implicit_t* cnav = (cnav_implicit_t*)context;

  // Take a step.
  double t_actual;
  int status = CVode(cnav->cvode, t + dt, cnav->U, &t_actual, CV_NORMAL);
  ASSERT(status != CV_MEM_NULL);
  ASSERT(status != CV_NO_MALLOC);
  ASSERT(status != CV_ILL_INPUT);
  ASSERT(status != CV_LINIT_FAIL);
  ASSERT(status != CV_LSOLVE_FAIL);
  ASSERT(status != CV_RHSFUNC_FAIL);

  if ((status != CV_SUCCESS) && (status != CV_TSTOP_RETURN) && 
      (status != CV_ROOT_RETURN) && (status != CV_FIRST_RHSFUNC_FAIL))
  {
    switch(status)
    {
      case CV_TOO_CLOSE:
        polymec_error("dt (%g) is too small.", dt);
        break;
      case CV_TOO_MUCH_WORK:
        polymec_error("Advance took too many internal steps.");
        break;
      case CV_TOO_MUCH_ACC:
        polymec_error("The integrator could not achieve the desired accuracy.");
        break;
      case CV_ERR_FAILURE:
        polymec_error("Too many failures during one internal step.");
        break;
      case CV_CONV_FAILURE:
        polymec_error("Too many convergence failures.");
        break;
      case CV_REPTD_RHSFUNC_ERR:
        polymec_error("Too many errors in the right hand side function.");
        break;
      case CV_UNREC_RHSFUNC_ERR:
        polymec_error("Integrator failed to recover from a right hand side error.");
        break;
      case CV_RTFUNC_FAIL:
        polymec_error("Integrator could not find a root.");
        break;
      default:
        polymec_error("The integrator failed.");
        break;
    }
  }

  // If t_actual > t + dt, interpolate backward.
  // FIXME
}

static void cnav_reconnect(void* context, mesh_t* new_mesh)
{
}

static void cnav_init(void* context, double t)
{
  cnav_implicit_t* cnav = (cnav_implicit_t*)context;

  // Make sure our solution vector is allocated.
  if (cnav->U != NULL)
  {
    N_Vector_free(&cnav->U);
    CVodeFree(&cnav->cvode);
  }
  cnav->cvode = CVodeCreate(CV_BDF, CV_NEWTON);
  CVodeSetUserData(cnav->cvode, cnav);
  cnav->U = N_Vector_new(cnav->comm, cnav->mesh->num_cells);

  // Initialize the solution.
  double* U = N_Vector_data(cnav->U);
  for (int c = 0; c < num_cells; ++c)
    st_func_eval(cnav->initial_cond, &cnav->mesh->cells[c].center, t, &U[c]);
  CVodeInit(cnav->cvode, compute_F, t, cnav->U);

  // We try GMRES with the preconditioner applied on the left.
  CVSpgmr(cnav->cvode, PREC_LEFT, 0);

  // Set the Jacobian-times-vector function.
  //CVSpilsSetJacTimesVecFn(cnav->cvode, JoV);

  // Set modified Gram-Schmidt orthogonalization 
  //CVSpilsSetGSType(cnav->cvode, MODIFIED_GS);

  // Set the preconditioner solve and setup functions.
  //CVSpilsSetPreconditioner(cvode_mem, Precond, PSolve);
}

static void cnav_save(void* context, io_interface_t* io, double t, int step)
{
  cnav_implicit_t* cnav = (cnav_implicit_t*)context;
  model_save(cnav->model);
}

static void cnav_dtor(void* ctx)
{
  cnav_implicit_t* cnav = (cnav_implicit_t*)context;
  if (cnav->U != NULL)
    N_Vector_free(&cnav->U);
  if (cnav->cvode != NULL)
    CVodeFree(&cnav->cvode);
  free(cnav);
}

model_t* cnav_implicit_model_new(int order,
                                 options_t* options)
{
  ASSERT(order >= 1);
  ASSERT(order <= 5);

  model_vtable vtable = { .read_input = cnav_read_input,
                          .init = cnav_init,
                          .max_dt = cnav_max_dt,
                          .advance = cnav_advance,
                          .save = cnav_save,
                          .dtor = cnav_dtor};
  cnav_implicit_t* cnav = malloc(sizeof(cnav_implicit_t));
  model_t* model = model_new("Implicit compressible Navier-Stokes", cnav, vtable, options);

  // Initialize the bookkeeping structures.
#if USE_MPI
  cnav->comm = MPI_COMM_WORLD;
#endif
  cnav->cvode = NULL;
  cnav->U = NULL;

  // Set up the saver.
  io_interface_t* saver = silo_io_new(cnav->comm, 0, false);
  model_set_saver(model, saver);

  return model;
}

model_t* create_cnav_implicit(int order,
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
  model_t* model = cnav_implicit_model_new(order, options);
  cnav_implicit_t* cnav = (cnav_implicit_t*)model_context(model);
  cnav->mesh = mesh;
  // FIXME
  return model;
}

#ifdef __cplusplus
}
#endif

