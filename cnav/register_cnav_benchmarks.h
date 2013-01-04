#ifndef POLYMEC_REGISTER_CNAV_BENCHMARKS_H
#define POLYMEC_REGISTER_CNAV_BENCHMARKS_H

#include "core/model.h"

#ifdef __cplusplus
extern "C" {
#endif

// Registers all the cnav benchmarks with the model.
void register_cnav_benchmarks(model_t* model);

#ifdef __cplusplus
}
#endif

#endif
