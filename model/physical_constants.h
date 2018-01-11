// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_PHYSICAL_CONSTANTS_H
#define POLYMEC_PHYSICAL_CONSTANTS_H

#include "core/polymec.h"

// This subsystem allows applications to define systems of units, and 
// associate with each system a set of physical constants. Automatic 
// conversions are provided between any two system of units for all common
// physical constants.
//
// Names of unit systems and constants are non-negative integers. They can 
// also be enumerated types (which can be converted to integers in C).
//
// Support is also offered for a globally-adopted, or "default" system of 
// units. By default, there is no valid system of units.

// Define a new system of units using the given label.
void define_units_system(int units_system_name);

// Set the global system of units to the give (defined) one. This results 
// in a run-time error if the given system of units has not been defined.
void use_units_system(int units_system_name);

// Defines a physical constant within a given system of units.
void set_physical_constant(int physical_constant_name, real_t value, int units_system_name);

// Retrieves a physical constant using the global system of units.
real_t physical_constant(int physical_constant_name);

// Retrieves a physical constant using the given system of units.
real_t physical_constant_in_units(int physical_constant_name, int units_system_name);

// Returns a conversion factor for a physical constant in moving from 
// units_system_1 to units_system_2.
real_t physical_constant_conversion_factor(int physical_constant_name, 
                                           int units_system_1,
                                           int units_system_2);

#endif
