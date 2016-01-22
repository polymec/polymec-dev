// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/unordered_map.h"
#include "model/physical_constants.h"

// Systems of units.
static int_ptr_unordered_map_t* units_systems = NULL;
static int_real_unordered_map_t* current_units = NULL;

void define_units_system(int units_system_name)
{
  if (units_systems == NULL)
    units_systems = int_ptr_unordered_map_new();
  if (int_ptr_unordered_map_contains(units_systems, units_system_name))
    polymec_error("define_units_system: units system %d is already defined.", units_system_name);
  int_ptr_unordered_map_insert_with_v_dtor(units_systems, 
                                           units_system_name, 
                                           int_real_unordered_map_new(),
                                           DTOR(int_real_unordered_map_free));
}

void use_units_system(int units_system_name)
{
  if (units_systems == NULL)
    polymec_error("use_units_system: No valid units system is defined!");
  int_real_unordered_map_t** ptr = (int_real_unordered_map_t**)int_ptr_unordered_map_get(units_systems, units_system_name);
  if (ptr == NULL)
    polymec_error("use_units_system: units system %d is not defined!", units_system_name);
  current_units = *ptr;
}

void set_physical_constant(int physical_constant_name, real_t value, int units_system_name)
{
  if (units_systems == NULL)
    polymec_error("set_physical_constant: No valid units system is defined!");
  int_real_unordered_map_t** ptr = (int_real_unordered_map_t**)int_ptr_unordered_map_get(units_systems, units_system_name);
  if (ptr == NULL)
    polymec_error("set_physical_constant: units system %d is not defined!", units_system_name);
  int_real_unordered_map_t* units = *ptr;
  int_real_unordered_map_insert(units, physical_constant_name, value);
}

real_t physical_constant(int physical_constant_name)
{
  ASSERT(current_units != NULL);
  real_t* val = int_real_unordered_map_get(current_units, physical_constant_name);
  ASSERT(val != NULL);
  return *val;
}

real_t physical_constant_in_units(int physical_constant_name, int units_system_name)
{
  ASSERT(units_systems != NULL);
  int_real_unordered_map_t** ptr = (int_real_unordered_map_t**)int_ptr_unordered_map_get(units_systems, units_system_name);
  ASSERT(ptr != NULL);
  int_real_unordered_map_t* units = *ptr;
  real_t* val = int_real_unordered_map_get(units, physical_constant_name);
  ASSERT(val != NULL);
  return *val;
}

real_t physical_constant_conversion_factor(int physical_constant_name, 
                                           int units_system_1,
                                           int units_system_2)
{
  ASSERT(units_systems != NULL);

  int_real_unordered_map_t** ptr1 = (int_real_unordered_map_t**)int_ptr_unordered_map_get(units_systems, units_system_1);
  ASSERT(ptr1 != NULL);
  int_real_unordered_map_t* units1 = *ptr1;
  real_t* val1 = int_real_unordered_map_get(units1, physical_constant_name);
  ASSERT(val1 != NULL);
  ASSERT(*val1 != 0.0);

  int_real_unordered_map_t** ptr2 = (int_real_unordered_map_t**)int_ptr_unordered_map_get(units_systems, units_system_2);
  ASSERT(ptr2 != NULL);
  int_real_unordered_map_t* units2 = *ptr2;
  real_t* val2 = int_real_unordered_map_get(units2, physical_constant_name);
  ASSERT(val2 != NULL);

  return (*val2)/(*val1);
}

