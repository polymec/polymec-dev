#ifndef POLYMEC_ADVECT_MODEL_H
#define POLYMEC_ADVECT_MODEL_H

#include "core/model.h"

// Creates an advection/diffusion/reaction model using the given options.
model_t* advect_model_new(options_t* options);

// This factory method creates a new advection model object that is ready 
// to run a problem defined by the given parameters.
model_t* create_advect(mesh_t* mesh,
                       st_func_t* velocity, 
                       st_func_t* diffusivity, 
                       st_func_t* source, 
                       st_func_t* initial_cond, 
                       string_ptr_unordered_map_t* bcs, 
                       st_func_t* solution,
                       options_t* options);


#endif

