#ifndef POLYMEC_SUNDIALS_HELPERS_H
#define POLYMEC_SUNDIALS_HELPERS_H

#if USE_MPI
#include "nvector/nvector_parallel.h"
#else
#include "nvector/nvector_serial.h"
#endif

// These macros and functions help with creating and manipulating serial 
// and parallel N_Vector objects.
#if USE_MPI
#define NV_DATA(v) NV_DATA_P(v)
#define NV_LOCLENGTH(v) NV_LOCLENGTH_P(v)
#define NV_GLOBLENGTH(v) NV_GLOBLENGTH_P(v)
#define NV_Ith(v) NV_Ith_P(v)
#else
#define NV_DATA(v) NV_DATA_S(v)
#define NV_LOCLENGTH(v) NV_LENGTH_S(v)
#define NV_GLOBLENGTH(v) NV_LENGTH_S(v)
#define NV_Ith(v) NV_Ith_S(v)
#endif

#endif


