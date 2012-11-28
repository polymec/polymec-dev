#ifndef POLYMEC_POISSON_H
#define POLYMEC_POISSON_H

#include "core/model.h"

#ifdef __cplusplus
extern "C" {
#endif

// Creates a Poisson model using the given options.
model_t* poisson_model_new(options_t* options);

#ifdef __cplusplus
}
#endif

#endif

