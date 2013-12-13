// Copyright (c) 2012-2013, Jeffrey N. Johnson
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this 
// list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice, 
// this list of conditions and the following disclaimer in the documentation 
// and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#ifndef POLYMEC_SUNDIALS_HELPERS_H
#define POLYMEC_SUNDIALS_HELPERS_H

#include "mpi.h"

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
#define NV_Ith(v,i) NV_Ith_P(v,i)
#else
#define NV_DATA(v) NV_DATA_S(v)
#define NV_LOCLENGTH(v) NV_LENGTH_S(v)
#define NV_GLOBLENGTH(v) NV_LENGTH_S(v)
#define NV_Ith(v,i) NV_Ith_S(v,i)
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

static inline N_Vector N_VNewEmpty(MPI_Comm comm, int local_len)
{
#if HAVE_MPI
  int global_len;
  MPI_Allreduce(&local_len, &global_len, 1, MPI_INT, MPI_SUM, comm);
  return N_VNewEmpty_Parallel(comm, local_len, global_len);
#else
  return N_VNewEmpty_Serial(local_len);
#endif
}

static inline N_Vector N_VMake(MPI_Comm comm, int local_len, double* v_data)
{
#if HAVE_MPI
  int global_len;
  MPI_Allreduce(&local_len, &global_len, 1, MPI_INT, MPI_SUM, comm);
  return N_VMake_Parallel(comm, local_len, global_len, v_data);
#else
  return N_VMake_Serial(local_len, v_data);
#endif
}

static inline void N_VPrint(N_Vector v)
{
#if HAVE_MPI
  N_VPrint_Parallel(v);
#else
  N_VPrint_Serial(v);
#endif
}

#endif


