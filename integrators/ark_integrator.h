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

#ifndef POLYMEC_ARK_INTEGRATOR_H
#define POLYMEC_ARK_INTEGRATOR_H

#include "integrators/integrator.h"

// Types of Additive Runge-Kutta integrators.
typedef enum
{
  ASIRK_1A,    // Semi-implicit 1st order with full (non)linear solve
  ASIRK_1B,    // Semi-implicit 1st order with implicit linearization 
  ASIRK_1C     // Semi-implicit 1st order with implicit linearization (variant)
} ark_integrator_type_t;

// Creates a semi-implicit Additive Runge-Kutta integrator of the given order,
// using the given explicit integrator to integrate non-stiff terms, and the 
// given implicit integrator to integrate stiff and diffusive terms.
integrator_t* ark_integrator_new(ark_integrator_type_t type,
                                 integrator_t* explicit_integrator,
                                 integrator_t* implicit_integrator);

#endif

