#ifndef POLYMEC_POISSON_H
#define POLYMEC_POISSON_H

#include "core/model.h"

// Creates a Poisson model using the given options.
model_t* poisson_model_new(options_t* options);

// This factory method creates a new Poisson model object that is ready 
// to run a problem defined by the given parameters.
model_t* create_poisson(mesh_t* mesh,
                        st_func_t* rhs,
                        str_ptr_unordered_map_t* bcs, 
                        st_func_t* solution,
                        options_t* options);

#endif

