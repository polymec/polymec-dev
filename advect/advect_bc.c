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

#include <gc/gc.h>
#include "advect/advect_bc.h"
#include "core/st_func.h"
#include "core/periodic_bc.h"

static void advect_bc_free(void* bc, void* dummy)
{
  if (pointer_is_periodic_bc(bc)) return;
  advect_bc_t* pbc = (advect_bc_t*)bc;
  pbc->F = NULL;
}

advect_bc_t* advect_bc_new(double alpha, double beta, st_func_t* F)
{
  ASSERT(F != NULL);
  advect_bc_t* bc = GC_MALLOC(sizeof(advect_bc_t));
  bc->alpha = alpha;
  bc->beta = beta;
  bc->F = F;
  GC_register_finalizer(bc, advect_bc_free, bc, NULL, NULL);
  return bc;
}

