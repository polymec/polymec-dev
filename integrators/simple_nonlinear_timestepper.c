#include "integrators/simple_nonlinear_timestepper.h"

typedef struct
{
  double initial_step_size;
  double max_step_size;
  double reduction_factor;
  int max_iterations;
  double increase_factor;
} simple_stepper_t;

static double simple_compute_dt(void* context,
                                double last_unsuccessful_dt,
                                int num_failures,
                                double* dt_history,
                                double* error_history,
                                int* iteration_history,
                                int history_length,
                                char* explanation)
{
  simple_stepper_t* stepper = context;
  double dt = dt_history[0];
  if (history_length == 0)
  {
    sprintf(explanation, "Initial step size is %g", stepper->initial_step_size);
    dt = stepper->initial_step_size;
  }
  else if (num_failures >= stepper->max_iterations)
  {
    sprintf(explanation, "Number of failures (%d) exceeded threshold.", num_failures);
    dt = stepper->reduction_factor * last_unsuccessful_dt;
  }
  else if (num_failures == 0)
  {
    sprintf(explanation, "Converged first time.");
    dt = stepper->increase_factor * dt_history[0];
  }
  else
    sprintf(explanation, "No change to step size.");

  if (dt > stepper->max_step_size)
  {
    sprintf(explanation, "Maximum allowable step size.");
    dt = stepper->max_step_size;
  }
  return dt;
}

static bool simple_recompute_J(void* context, 
                               double* dt_history,
                               double* error_history,
                               int* iteration_history,
                               int history_length)
{
  return true;
}

static void simple_dtor(void* context)
{
  simple_stepper_t* stepper = context;
  free(stepper);
}

nonlinear_timestepper_t* simple_nonlinear_timestepper_new(double initial_step_size,
                                                          double max_step_size,
                                                          double reduction_factor,
                                                          int max_iterations,
                                                          double increase_factor)
{
  ASSERT(initial_step_size > 0.0);
  ASSERT(max_step_size >= initial_step_size);
  ASSERT(reduction_factor > 0.0);
  ASSERT(reduction_factor < 1.0);
  ASSERT(max_iterations > 1);
  ASSERT(increase_factor > 1.0);

  char name[1024];
  snprintf(name, 1024, "Simple timestepper (max step size = %g, reduction factor = %g,\n"
                       "max iterations = %d, increase factor = %g)", 
           max_step_size, reduction_factor, max_iterations, increase_factor);
  nonlinear_timestepper_vtable vtable = {.compute_dt = simple_compute_dt,
                                         .recompute_J = simple_recompute_J,
                                         .dtor = simple_dtor};
  simple_stepper_t* context = malloc(sizeof(simple_stepper_t));
  context->initial_step_size = initial_step_size;
  context->max_step_size = max_step_size;
  context->reduction_factor = reduction_factor;
  context->max_iterations = max_iterations;
  context->increase_factor = increase_factor;
  return nonlinear_timestepper_new(name, context, 1, vtable);
}

