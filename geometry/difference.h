#ifndef POLYMEC_DIFFERENCE_H
#define POLYMEC_DIFFERENCE_H

#include "core/sp_func.h"

#ifdef __cplusplus
extern "C" {
#endif

// This signed distance function represents the difference of two 
// given surfaces represented by signed distance functions.
sp_func_t* difference_new(sp_func_t* surface1, sp_func_t* surface2);

#ifdef __cplusplus
}
#endif

#endif

