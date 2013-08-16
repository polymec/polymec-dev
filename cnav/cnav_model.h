#ifndef POLYMEC_CNAV_MODEL_H
#define POLYMEC_CNAV_MODEL_H

#include "core/model.h"
#include "integrators/integrator.h"
#include "cnav/cnav_eos.h"

// Types of time integrators available for this model.
typedef enum
{
  CNAV_SEMI_IMPLICIT,  // Explicit advection, implicit diffusion/reaction
  CNAV_IMPLICIT        // Fully-coupled implicit integration
} cnav_time_integrator_t;

// Creates a compressible Navier-Stokes model using the given options.
model_t* cnav_model_new(options_t* options);

// This factory method creates a new compressible Navier-Stokes model object 
// that is ready to run a problem defined by the given parameters.
model_t* create_cnav(cnav_time_integrator_t integrator,
                     int order,
                     mesh_t* mesh,
                     cnav_eos_t* equation_of_state,
                     st_func_t* source,
                     st_func_t* initial_cond, 
                     str_ptr_unordered_map_t* bcs, 
                     st_func_t* solution,
                     options_t* options);


#endif

