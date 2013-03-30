#ifndef POLYMEC_SUNDIALS_HELPERS_H
#define POLYMEC_SUNDIALS_HELPERS_H

#if HAVE_MPI
#include "nvector/nvector_parallel.h"
#else
#include "nvector/nvector_serial.h"
#endif
#include "sundials/sundials_iterative.h"

// These macros and functions help with creating and manipulating serial 
// and parallel N_Vector objects.
#if HAVE_MPI
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

static inline N_Vector N_VNew(MPI_Comm comm, int local_len)
{
#if HAVE_MPI
  int global_len;
  MPI_Allreduce(&local_len, &global_len, 1, MPI_INT, MPI_SUM, comm);
  return N_VNew_Parallel(comm, local_len, global_len);
#else
  return N_VNew_Serial(local_len);
#endif
}

#endif


