#ifndef ARBI_SIMULATION_H
#define ARBI_SIMULATION_H

#include "arbi.h"
#include "options.h"
#include "model.h"
#include "soln_vector.h"

#ifdef __cplusplus
extern "C" {
#endif

// A simulation represents a numerical problem to be solved.
typedef struct simulation_t simulation_t;

// Construct a new simulation with the given model and options.
simulation_t* simulation_new(model_t* model, options_t* options);

// Destroy the simulation.
void simulation_free(simulation_t* sim);

// Returns the (simulation) time at which the simulation commences.
double simulation_start_time(simulation_t* sim);

// Returns the maximum number of steps for which the simulation will run.
int simulation_num_steps(simulation_t* sim);

// Returns the maximum time to which the simulation will run.
double simulation_max_time(simulation_t* sim);

// Returns the model bound to the simulation.
model_t* simulation_model(simulation_t* sim);

// Initializes the given simulation and solution vector at the given time.
void simulation_init(simulation_t* sim, double t);

// Advances the solution in the simulation by a single time step.
void simulation_advance(simulation_t* sim);

// Returns the current simulation time.
double simulation_time(simulation_t* sim);

// Returns the current simulation step.
int simulation_step(simulation_t* sim);

// Runs the current simulation to completion.
void simulation_run(simulation_t* sim);

#ifdef __cplusplus
}
#endif

#endif

