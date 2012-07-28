#ifndef ARBI_SIMULATION_H
#define ARBI_SIMULATION_H

#include "arbi.h"
#include "options.h"
#include "soln_vector.h"

#ifdef __cplusplus
extern "C" {
#endif

// A simulation represents a numerical problem to be solved.
typedef struct simulation_t simulation_t;

// Construct a new simulation from the data in an input file, along with 
// any command-line arguments.
simulation_t* simulation_from_file(FILE* file, options_t* opts);

// Destroy the simulation.
void simulation_free(simulation_t* sim);

// Create a solution vector to store the solution for this simulation.
soln_vector_t* simulation_create_vector(simulation_t* sim);

// Returns the (simulation) time at which the simulation commences.
double simulation_start_time(simulation_t* sim);

// Returns the maximum number of steps for which the simulation will run.
int simulation_num_steps(simulation_t* sim);

// Returns the maximum time to which the simulation will run.
double simulation_max_time(simulation_t* sim);

// Initializes the given simulation and solution vector at the given time.
int simulation_init(simulation_t* sim, soln_vector_t* solution, double t);

// Invokes all callbacks (periodic work, etc) to be called by the simulation
// at the given simulation time and/or step.
void simulation_invoke_callbacks(simulation_t* sim, soln_vector_t* solution, double t, int step);

// Advances the solution in the simulation by a single time step.
int simulation_step(simulation_t* sim, soln_vector_t* solution, double* t);

#ifdef __cplusplus
}
#endif

#endif

