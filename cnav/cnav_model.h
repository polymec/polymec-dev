#ifndef POLYMEC_CNAV_MODEL_H
#define POLYMEC_CNAV_MODEL_H

#include "core/model.h"

#ifdef __cplusplus
extern "C" {
#endif

// Creates an cnavion/diffusion/reaction model using the given options.
model_t* cnav_model_new(options_t* options);

// This factory method creates a new cnavion model object that is ready 
// to run a problem defined by the given parameters.
model_t* create_cnav(mesh_t* mesh,
                     st_func_t* initial_cond, 
                     str_ptr_unordered_map_t* bcs, 
                     options_t* options);


#ifdef __cplusplus
}
#endif

#endif

