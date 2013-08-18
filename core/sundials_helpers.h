// Copyright 2012-2013 Jeffrey Johnson.
// 
// This file is part of Polymec, and is licensed under the Apache License, 
// Version 2.0 (the "License"); you may not use this file except in 
// compliance with the License. You may may find the text of the license in 
// the LICENSE file at the top-level source directory, or obtain a copy of 
// it at
// 
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

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


