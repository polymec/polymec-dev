#ifndef POLYMEC_CONSTANT_ST_FUNC_H
#define POLYMEC_CONSTANT_ST_FUNC_H

#include "core/st_func.h"

#ifdef __cplusplus
extern "C" {
#endif

// Construct a constant space-time function given a number of components 
// and their values.
st_func_t* constant_st_func_new(int num_comp, double comp[]);

// Free of charge, we toss in the sp_func version.
sp_func_t* constant_sp_func_new(int num_comp, double comp[]);

#ifdef __cplusplus
}
#endif

#endif

