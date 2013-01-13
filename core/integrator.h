#ifndef POLYMEC_INTEGRATOR_H
#define POLYMEC_INTEGRATOR_H

#include "core/polymec.h"

#ifdef __cplusplus
extern "C" {
#endif

// This class provides an abstract interface for integrating systems of 
// ordinary differential equations. 
typedef struct integrator_t integrator_t;

// A function for initializing a given solution at time t.
typedef void (*integrator_init_func)(void*, double, double*, int);

// A function for integrating a given solution from time t1 to t2.
// The solution is integrated in place and has dimension N.
typedef void (*integrator_step_func)(void*, double, double, double*, int);

// A destructor function for the context object (if any).
typedef void (*integrator_dtor)(void*);

// This virtual table must be implemented by any integrator.
typedef struct 
{
  integrator_step_func step;
  integrator_init_func init;
  integrator_dtor      dtor;
} integrator_vtable;

// Types of integrators (for algorithmic classification).
typedef enum
{
  INTEGRATOR_EXPLICIT,
  INTEGRATOR_IMPLICIT,
  INTEGRATOR_IMEX
} integrator_type_t;

// Creates an integrator solver with the given name, context, and virtual table.
// Also given are some metadata, such as the order of the method, whether 
// it is explicit, implicit, or hybrid.
integrator_t* integrator_new(const char* name, 
                             void* context,
                             integrator_vtable vtable,
                             int order,
                             integrator_type_t type);

// Frees an integrator.
void integrator_free(integrator_t* integrator);

// Returns the name of the integrator (internally stored).
char* integrator_name(integrator_t* integrator);

// Returns the context object for this integrator.
void* integrator_context(integrator_t* integrator);

// Returns the order of the integration method.
int integrator_order(integrator_t* integrator);

// Returns the type of the integrator (explicit, implicit, IMEX).
integrator_type_t integrator_type(integrator_t* integrator);

// Initializes the integrator with a solution at the given time 
// (and a given number of unknowns).
void integrator_init(integrator_t* integrator, double t, double* solution, int N);

// Integrates the given solution from time t1 to t2.
// The solution is integrated in place.
void integrator_step(integrator_t* integrator, double t1, double t2, 
                     double* solution);

#ifdef __cplusplus
}
#endif

#endif

