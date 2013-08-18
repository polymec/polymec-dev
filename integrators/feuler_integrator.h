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

#ifndef POLYMEC_FEULER_INTEGRATOR_H
#define POLYMEC_FEULER_INTEGRATOR_H

#include "integrators/integrator.h"

// A function signature for explicitly computing the derivative of a solution 
// at the given time. Arguments are:
// 1. A context object
// 2. The time t
// 3. The solution u 
// 4. Storage for the derivative of the solution, du/dt.
typedef void (*feuler_compute_deriv)(void*, double, double*, double*);

// Creates an integrator that uses the 1st-order explicit forward Euler method,
// using the given function to compute the derivative of the solution at t1.
integrator_t* feuler_integrator_new(void* context, 
                                    feuler_compute_deriv compute_deriv,
                                    integrator_dtor dtor);

#endif

