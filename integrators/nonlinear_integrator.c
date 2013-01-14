#include "integrators/nonlinear_integrator.h"
#include "cvode/cvode_spgmr.h"
#include "cvode/cvode_spbcgs.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
  MPI_Comm comm;
  void* context;
  void* cvode;
  int type, order;
  N_Vector u;
  integrator_compute_F_func compute_F;
  integrator_compute_Jv_func compute_Jv;
  nonlinear_integrator_solver_type_t solver_type;
  int gram_schmidt, max_kdim;
  int precond_type;
  nonlinear_integrator_precond_setup_func precond_setup;
  nonlinear_integrator_precond_solve_func precond_solve;
  integrator_dtor dtor;

  bool stepped;
} nonlinear_integrator_t;

static void nl_reset(nonlinear_integrator_t* nli)
{
  if (nli->cvode != NULL)
    N_VDestroy(nli->u);
  nli->u = NULL;
  nli->stepped = false;
}

static void nl_init(void* context, int N)
{
  nonlinear_integrator_t* nli = context;
  nl_reset(nli);
  nli->u = N_VNew(nli->comm, N);
}

static void nl_step(void* context, double t1, double t2, double* solution, int N)
{
  ASSERT(t2 > t1);
  double dt = t2 - t1;

  nonlinear_integrator_t* nli = context;
  ASSERT(nli->cvode != NULL);
  if (!nli->stepped)
  {
    CVodeInit(nli->cvode, nli->compute_F, t1, nli->u);
    nli->stepped = true;
  }

  // Take a step.
  double t_actual;
  int status = CVode(nli->cvode, t1 + dt, nli->u, &t_actual, CV_NORMAL);
  ASSERT(status != CV_MEM_NULL);
  ASSERT(status != CV_NO_MALLOC);
  ASSERT(status != CV_ILL_INPUT);
  ASSERT(status != CV_LINIT_FAIL);
  ASSERT(status != CV_LSOLVE_FAIL);
  ASSERT(status != CV_RHSFUNC_FAIL);

  if ((status != CV_SUCCESS) && (status != CV_TSTOP_RETURN) && 
      (status != CV_ROOT_RETURN)) // && (status != CV_FIRST_RHSFUNC_FAIL))
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

  // FIXME: Copy solution to array!
}

static void nl_dtor(void* context)
{
  nonlinear_integrator_t* nli = context;
  nl_reset(nli);
  if (nli->cvode != NULL)
    CVodeFree(nli->cvode);
  if ((nli->context != NULL) && (nli->dtor != NULL))
    nli->dtor(nli->context);
  free(nli);
}


integrator_t* nonlinear_integrator_new(MPI_Comm comm,
                                       void* context, 
                                       int cvode_type,
                                       int order,
                                       integrator_compute_F_func compute_F,
                                       integrator_compute_Jv_func compute_Jv,
                                       nonlinear_integrator_solver_type_t solver_type,
                                       int max_kdim,
                                       int gram_schmidt,
                                       int precond_type,
                                       nonlinear_integrator_precond_setup_func precond_setup,
                                       nonlinear_integrator_precond_solve_func precond_solve,
                                       integrator_dtor dtor)
{
  ASSERT((cvode_type == CV_ADAMS) || (cvode_type == CV_BDF));
  ASSERT(order <= 12);
  ASSERT((order <= 5) || (cvode_type == CV_ADAMS));
  ASSERT(max_kdim > 1);
  ASSERT((gram_schmidt == CLASSICAL_GS) || (gram_schmidt == MODIFIED_GS));
  ASSERT((precond_type == PREC_NONE) || (precond_type == PREC_LEFT) ||
         (precond_type == PREC_RIGHT) || (precond_type == PREC_BOTH));
  nonlinear_integrator_t* nli = malloc(sizeof(nonlinear_integrator_t));
  nli->comm = comm;
  nli->context = context;
  nli->type = cvode_type;
  nli->order = order;
  nli->compute_F = compute_F;
  nli->compute_Jv = compute_Jv;
  nli->solver_type = solver_type;
  nli->max_kdim = max_kdim;
  nli->gram_schmidt = gram_schmidt;
  nli->precond_type = precond_type;
  nli->precond_setup = precond_setup;
  nli->precond_solve = precond_solve;
  nli->dtor = dtor;
  nl_reset(nli);
  if (nli->type == CV_ADAMS)
    nli->cvode = CVodeCreate(CV_ADAMS, CV_FUNCTIONAL);
  else
    nli->cvode = CVodeCreate(CV_BDF, CV_NEWTON);
  CVodeSetUserData(nli->cvode, nli);
  if (nli->solver_type == INTEGRATOR_SOLVER_GMRES)
  {
    CVSpgmr(nli->cvode, nli->precond_type, nli->max_kdim);
    CVSpilsSetGSType(nli->cvode, nli->gram_schmidt);
  }
  else
  {
    CVSpbcg(nli->cvode, nli->precond_type, nli->max_kdim);
    CVSpilsSetGSType(nli->cvode, nli->gram_schmidt);
  }
  CVSpilsSetJacTimesVecFn(nli->cvode, nli->compute_Jv);
  CVSpilsSetPreconditioner(nli->cvode, nli->precond_setup, nli->precond_solve);

  integrator_vtable vtable = { .init = nl_init, .step = nl_step, .dtor = nl_dtor };
  char name[1024];
  if (cvode_type == CV_ADAMS)
    sprintf(name, "Adams (non-stiff, order %d)", order);
  else
    sprintf(name, "BDF (stiff, order %d)", order);
  return integrator_new(name, nli, vtable, order, INTEGRATOR_IMPLICIT);
}

#ifdef __cplusplus
}
#endif

