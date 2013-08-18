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

#ifndef POLYMEC_CNAV_MODEL_H
#define POLYMEC_CNAV_MODEL_H

#include "core/model.h"
#include "integrators/integrator.h"
#include "cnav/cnav_eos.h"

// Types of time integrators available for this model.
typedef enum
{
  CNAV_SEMI_IMPLICIT,  // Explicit advection, implicit diffusion/reaction
  CNAV_IMPLICIT        // Fully-coupled implicit integration
} cnav_time_integrator_t;

// Creates a compressible Navier-Stokes model using the given options.
model_t* cnav_model_new(options_t* options);

// This factory method creates a new compressible Navier-Stokes model object 
// that is ready to run a problem defined by the given parameters.
model_t* create_cnav(cnav_time_integrator_t integrator,
                     int order,
                     mesh_t* mesh,
                     cnav_eos_t* equation_of_state,
                     st_func_t* source,
                     st_func_t* initial_cond, 
                     string_ptr_unordered_map_t* bcs, 
                     st_func_t* solution,
                     options_t* options);


#endif

