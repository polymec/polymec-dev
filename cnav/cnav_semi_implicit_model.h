#ifndef POLYMEC_CNAV_SEMI_IMPLICIT_MODEL_H
#define POLYMEC_CNAV_SEMI_IMPLICIT_MODEL_H

#include "core/model.h"

#ifdef __cplusplus
extern "C" {
#endif

// Creates a semi-implicit compressible Navier-Stokes model using the 
// given options.
model_t* cnav_semi_implicit_model_new(options_t* options);

// This factory method creates a new cnavion model object that is ready 
// to run a problem defined by the given parameters.
model_t* create_cnav_semi_implicit(mesh_t* mesh,
                                   cnav_eos_t* equation_of_state,
//                                  reaction_network_t* reactions,
                                   st_func_t* source, 
                                   st_func_t* initial_cond, 
                                   str_ptr_unordered_map_t* bcs, 
                                   st_func_t* solution,
                                   options_t* options);


#ifdef __cplusplus
}
#endif

#endif

