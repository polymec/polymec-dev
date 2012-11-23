#ifndef ARBI_ADVECT_MODEL_H
#define ARBI_ADVECT_MODEL_H

#include "core/model.h"

#ifdef __cplusplus
extern "C" {
#endif

// Creates an advection/diffusion/reaction model using the given options.
model_t* advect_model_new(options_t* options);

#ifdef __cplusplus
}
#endif

#endif

