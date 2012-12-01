#ifndef POLYMEC_REGISTER_ADVECT_BENCHMARKS_H
#define POLYMEC_REGISTER_ADVECT_BENCHMARKS_H

#include "core/model.h"

#ifdef __cplusplus
extern "C" {
#endif

// Registers all the advect benchmarks with the model.
void register_advect_benchmarks(model_t* model);

#ifdef __cplusplus
}
#endif

#endif
