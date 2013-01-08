#ifndef POLYMEC_CNAV_IMPLICIT_MODEL_H
#define POLYMEC_CNAV_IMPLICIT_MODEL_H

#include "core/unordered_map.h"
#include "core/st_func.h"
#include "core/model.h"
#include "cnav/cnav_eos.h"

#ifdef __cplusplus
extern "C" {
#endif

// Creates an implicit compressible Navier-Stokes model using the given 
// time integrator and options.
model_t* cnav_implicit_model_new(int order, options_t* options);

// This factory method creates a new implicit compressible Navier-Stokes 
// model object that is ready to run a problem defined by the given parameters.
model_t* create_cnav_implicit(int order,
                              mesh_t* mesh,
                              cnav_eos_t* equation_of_state,
                              st_func_t* source,
                              st_func_t* initial_cond, 
                              str_ptr_unordered_map_t* bcs, 
                              st_func_t* solution,
                              options_t* options);


#ifdef __cplusplus
}
#endif

#endif

