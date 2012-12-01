#ifndef POLYMEC_REGISTER_POISSON_BENCHMARKS_H
#define POLYMEC_REGISTER_POISSON_BENCHMARKS_H

#include "core/model.h"

#ifdef __cplusplus
extern "C" {
#endif

// Registers all the poisson benchmarks with the model.
void register_poisson_benchmarks(model_t* model);

#ifdef __cplusplus
}
#endif

#endif
