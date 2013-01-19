#ifndef POLYMEC_INTEGRATOR_H
#define POLYMEC_INTEGRATOR_H

#include "core/polymec.h"
#include "sundials/sundials_iterative.h"

#ifdef __cplusplus
extern "C" {
#endif

// The following types and functions are provided for easy interoperability
// with CVODE.

// A prototype for a function that computes a matrix-vector product for a
// linear system.
typedef int (*integrator_compute_Ax_func)(void *A_data, N_Vector x, N_Vector Ax);

// A prototype for a function that computes a function for a
// nonlinear system.
typedef int (*integrator_compute_F_func)(double t, N_Vector u, N_Vector u_dot, void* context);

// A prototype for a function that computes a matrix-vector product for a
// nonlinear system.
typedef int (*integrator_compute_Jv_func)(N_Vector v, N_Vector Jv, double t, N_Vector y, N_Vector fy, void* context, N_Vector tmp);

// A prototype for a function that computes right-hand-side vectors for 
// implicit integrators at a given time.
typedef void (*integrator_compute_rhs_func)(void*, double, N_Vector);

// A prototype for a function that applies discrete boundary conditions to 
// a linear system at a given time.
typedef void (*integrator_apply_bcs_func)(void*, double, N_Vector);

// A prototype for a function that solves the preconditioner system
// M * z = r. Here, precond_type = PRECOND_LEFT for left-preconditioned 
// systems and PRECOND_RIGHT for right-preconditioned systems.
typedef int (*integrator_precond_func)(void *P_data, N_Vector r, N_Vector z, int precond_type);


// This class provides an abstract interface for integrating systems of 
// ordinary differential equations. 
typedef struct integrator_t integrator_t;

// A function for initializing an integrator with a given dimension.
typedef void (*integrator_init_func)(void*, int);

// A function for integrating a given solution from time t1 to t2.
// The solution is integrated in place and has dimension N.
typedef void (*integrator_step_func)(void*, double, double, double*, int);

// A destructor function for the context object (if any).
typedef void (*integrator_dtor)(void*);

// This virtual table must be implemented by any integrator.
typedef struct 
{
  integrator_step_func       step;
  integrator_init_func       init;
  integrator_compute_Jv_func compute_Jv; // Optional
  integrator_dtor            dtor;
} integrator_vtable;

// Types of integrators (for algorithmic classification).
typedef enum
{
  INTEGRATOR_EXPLICIT,
  INTEGRATOR_IMPLICIT,
  INTEGRATOR_SEMI_IMPLICIT
} integrator_type_t;

// Creates an integrator with the given name, context, and virtual table.
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

// Initializes the integrator with a given number of unknowns.
void integrator_init(integrator_t* integrator, int N);

// Integrates the given solution from time t1 to t2.
// The solution is integrated in place.
void integrator_step(integrator_t* integrator, double t1, double t2, 
                     double* solution);

// Computes the product of the Jacobian matrix associated with this 
// (nonlinear implicit) integrator with the given vector. This is only 
// callable on integrators for which the compute_Jv entry is given in the 
// virtual table. Call only if you know what you're doing!
void integrator_compute_Jv(integrator_t* integrator, double t, N_Vector x, N_Vector F, N_Vector v, N_Vector Jv);

#ifdef __cplusplus
}
#endif

#endif

