#include <stdlib.h>
#include <float.h>
#include "simulation.h"

#ifdef __cplusplus
extern "C" {
#endif

struct simulation_t 
{
  double t0;     // Start time.
  double tmax;   // Max time.
  double t;      // Current time.
  int max_steps; // Max number of steps.
  int step;      // Current step.
  model_t* model;// Simulation model.
  char dt_expl[ARBI_MODEL_MAXDT_REASON_SIZE]; // Explanation of latest dt.
};

simulation_t* simulation_new(model_t* model, options_t* options)
{
  ASSERT(model != NULL);
  simulation_t* s = malloc(sizeof(simulation_t));
  s->t0 = s->t = 0.0;
  s->step = 0;
  s->tmax = FLT_MAX;
  s->max_steps = INT_MAX;
  s->model = model;
  return s;
}

void simulation_free(simulation_t* sim)
{
  // We don't own the model, so we don't delete it.
  free(sim);
}

double simulation_start_time(simulation_t* sim)
{
  return sim->t0;
}

int simulation_num_steps(simulation_t* sim)
{
  return sim->max_steps;
}

double simulation_max_time(simulation_t* sim)
{
  return sim->tmax;
}

double simulation_time(simulation_t* sim)
{
  return sim->t;
}

void simulation_init(simulation_t* sim, double t)
{
  sim->t = t;
  model_init(sim->model, t);
}

static void simulation_invoke_callbacks(simulation_t* sim)
{
}

int simulation_step(simulation_t* sim)
{
  return sim->step;
}

void simulation_advance(simulation_t* sim)
{
  ASSERT(sim->model != NULL);
  double dt = model_max_dt(sim->model, sim->t, sim->dt_expl);
  model_advance(sim->model, sim->t, dt);
  sim->t += dt;
}

model_t* simulation_model(simulation_t* sim)
{
  return sim->model;
}

void simulation_run(simulation_t* sim)
{
  while ((sim->t < sim->tmax) && (sim->step < sim->max_steps))
    simulation_advance(sim);
}
#ifdef __cplusplus
}
#endif

