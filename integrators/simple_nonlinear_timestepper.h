#ifndef POLYMEC_SIMPLE_NONLINEAR_TIMESTEPPER_H
#define POLYMEC_SIMPLE_NONLINEAR_TIMESTEPPER_H

#include "integrators/nonlinear_solver.h"

// Creates a simple timestepper that reduces the timestep by a reduction 
// factor upon failure to converge the given maximum number of iterations, 
// and increases it by an increase factor if the iteration converges the 
// first time. It signals to recompute the Jacobian after every successful
// nonlinear iteration.
nonlinear_timestepper_t* simple_nonlinear_timestepper_new(double initial_step_size,
                                                          double max_step_size,
                                                          double reduction_factor,
                                                          int max_iterations,
                                                          double increase_factor);
#endif

