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

#ifndef POLYMEC_SIMPLE_NONLINEAR_TIMESTEPPER_H
#define POLYMEC_SIMPLE_NONLINEAR_TIMESTEPPER_H

#include "integrators/nonlinear_solver.h"

// Creates a simple timestepper that reduces the timestep by a reduction 
// factor upon failure to converge the given maximum number of iterations, 
// and increases it by an increase factor if the iteration converges the 
// first time. It signals to recompute the Jacobian after every successful
// nonlinear iteration.
nonlinear_timestepper_t* simple_nonlinear_timestepper_new(double initial_step_size,
                                                          double max_step_size,
                                                          double reduction_factor,
                                                          int max_iterations,
                                                          double increase_factor);
#endif

