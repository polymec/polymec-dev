#include <stdlib.h>
#include "simulation.h"

#ifdef __cplusplus
extern "C" {
#endif

struct simulation_t 
{
  double t0;     // Start time.
  double tmax;   // Max time.
  int max_steps; // Max number of steps.
};

//------------------------------------------------------------------------
simulation_t* simulation_from_file(FILE* file, options_t* opts)
{
  simulation_t* s = malloc(sizeof(simulation_t));
  return s;
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
void simulation_free(simulation_t* sim)
{
  free(sim);
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
soln_vector_t* simulation_create_vector(simulation_t* sim)
{
  return NULL;
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
double simulation_start_time(simulation_t* sim)
{
  return sim->t0;
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
int simulation_num_steps(simulation_t* sim)
{
  return sim->max_steps;
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
double simulation_max_time(simulation_t* sim)
{
  return sim->tmax;
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
int simulation_init(simulation_t* sim, soln_vector_t* solution, double t)
{
  return CLIDE_SUCCESS;
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
void simulation_invoke_callbacks(simulation_t* sim, soln_vector_t* solution, double t, int step)
{
}
//------------------------------------------------------------------------

//------------------------------------------------------------------------
int simulation_step(simulation_t* sim, soln_vector_t* solution, double* t)
{
  return CLIDE_SUCCESS;
}
//------------------------------------------------------------------------

#ifdef __cplusplus
}
#endif

