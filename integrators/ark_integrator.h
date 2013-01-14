#ifndef POLYMEC_ARK_INTEGRATOR_H
#define POLYMEC_ARK_INTEGRATOR_H

#include "core/integrator.h"

#ifdef __cplusplus
extern "C" {
#endif

// Types of Additive Runge-Kutta integrators.
typedef enum
{
  ARK1_LINEAR, // Semi-implicit 1st order with linear implicit integrator
  ASIRK_1A,    // Semi-implicit 1st order with full nonlinear solve
  ASIRK_1B,    // Semi-implicit 1st order with implicit linearization 
  ASIRK_1C     // Semi-implicit 1st order with implicit linearization (variant)
} ark_integrator_type_t;

// Creates a semi-implicit Additive Runge-Kutta integrator of the given order,
// using the given explicit integrator to integrate non-stiff terms, and the 
// given implicit integrator to integrate stiff and diffusive terms.
integrator_t* ark_integrator_new(ark_integrator_type_t type,
                                 integrator_t* explicit_integrator,
                                 integrator_t* implicit_integrator);

#ifdef __cplusplus
}
#endif

#endif

