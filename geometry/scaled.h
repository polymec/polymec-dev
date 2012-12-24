#ifndef POLYMEC_SCALED_H
#define POLYMEC_SCALED_H

#include "core/sp_func.h"

#ifdef __cplusplus
extern "C" {
#endif

// This signed distance function takes another such function and scales it 
// by a given factor.
sp_func_t* scaled_new(sp_func_t* func, double scale_factor);

#ifdef __cplusplus
}
#endif

#endif

