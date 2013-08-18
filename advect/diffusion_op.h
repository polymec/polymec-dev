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

#ifndef POLYMEC_DIFFUSION_OP_H
#define POLYMEC_DIFFUSION_OP_H

#include "core/lin_op.h"
#include "core/st_func.h"

// Returns a 2nd-order, finite-volume, cell-centered diffusion operator.
lin_op_t* diffusion_op_new(mesh_t* mesh, st_func_t* diffusivity);

// Sets the time on the diffusion operator.
void diffusion_op_set_time(lin_op_t* op, double t);

#endif
