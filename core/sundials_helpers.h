// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_SUNDIALS_HELPERS_H
#define POLYMEC_SUNDIALS_HELPERS_H

#include "core/polymec.h"

#if POLYMEC_HAVE_MPI
#include "nvector/nvector_parallel.h"
#else
#include "nvector/nvector_serial.h"
#endif
#include "sundials/sundials_iterative.h"

// These macros and functions help with creating and manipulating serial
// and parallel N_Vector objects.
#if POLYMEC_HAVE_MPI
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

/// \collective Collective on comm.
static inline N_Vector N_VNew(MPI_Comm comm, int local_len)
{
#if POLYMEC_HAVE_MPI
  int global_len;
  MPI_Allreduce(&local_len, &global_len, 1, MPI_INT, MPI_SUM, comm);
  return N_VNew_Parallel(comm, local_len, global_len);
#else
  return N_VNew_Serial(local_len);
#endif
}

/// \collective Collective on comm.
static inline N_Vector N_VNewEmpty(MPI_Comm comm, int local_len)
{
#if POLYMEC_HAVE_MPI
  int global_len;
  MPI_Allreduce(&local_len, &global_len, 1, MPI_INT, MPI_SUM, comm);
  return N_VNewEmpty_Parallel(comm, local_len, global_len);
#else
  return N_VNewEmpty_Serial(local_len);
#endif
}

/// \collective Collective on comm.
static inline N_Vector N_VMake(MPI_Comm comm, int local_len, real_t* v_data)
{
#if POLYMEC_HAVE_MPI
  int global_len;
  MPI_Allreduce(&local_len, &global_len, 1, MPI_INT, MPI_SUM, comm);
  return N_VMake_Parallel(comm, local_len, global_len, v_data);
#else
  return N_VMake_Serial(local_len, v_data);
#endif
}

static inline void N_VPrint(N_Vector v)
{
#if POLYMEC_HAVE_MPI
  N_VPrint_Parallel(v);
#else
  N_VPrint_Serial(v);
#endif
}

#endif


