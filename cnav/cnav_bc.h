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

#ifndef POLYMEC_CNAV_BC_H
#define POLYMEC_CNAV_BC_H

#include <stdlib.h>
#include "core/st_func.h"

// Boundary condition structure for the compressible Navier-Stokes solver.
// Objects of this type are garbage-collected.
typedef struct
{
  double u_normal;       // Normal velocity (u dot n) at the boundary.
  vector_t u_tangential; // Tangential velocity at the boundary.

  double T_alpha, T_beta; // Coefficients for temperature (Robins) BCs.
} cnav_bc_t;

// Constructor for a cnav BC.
cnav_bc_t* cnav_bc_new();

#endif
